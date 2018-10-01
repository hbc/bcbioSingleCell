# FIXME Either make `sampleName` required or strip it from minimal examples.



# bcbioSingleCell ==============================================================
#' Read bcbio Single-Cell RNA-Seq Data
#'
#' @section Remote Data:
#'
#' When working in RStudio, we recommend connecting to the bcbio-nextgen run
#' directory as a remote connection over
#' [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @family S4 Generators
#' @author Michael Steinbaugh, Rory Kirchner
#' @export
#'
#' @inheritParams basejump::makeSummarizedExperiment
#' @inheritParams general
#' @param uploadDir `string`. Path to final upload directory. This path is set
#'   when running "`bcbio_nextgen -w template`".
#' @param sampleMetadataFile `string` or `NULL`. Sample barcode metadata file.
#'   Optional for runs with demultiplixed index barcodes (e.g. SureCell) and can
#'   be left `NULL`, but otherwise suggested for runs with multipliexed FASTQs
#'   containing multiple index barcodes (e.g. inDrops).
#' @param organism `string` or `NULL`. Organism name. Use the full Latin name
#'   (e.g. "Homo sapiens"), since this will be input downstream to AnnotationHub
#'   and ensembldb, unless `gffFile` is set. If left `NULL` (*not recommended*),
#'   the function call will skip loading gene-level annotations into
#'   [rowRanges()]. This can be useful for poorly annotation genomes or
#'   experiments involving multiple genomes.
#' @param ensemblRelease `scalar integer` or `NULL`. Ensembl release version. If
#'   `NULL`, defaults to current release, and does not typically need to be
#'   user-defined. Passed to AnnotationHub for `EnsDb` annotation matching,
#'   unless `gffFile` is set.
#' @param genomeBuild `string` or `NULL`. Ensembl genome build name (e.g.
#'   "GRCh38"). This will be passed to AnnotationHub for `EnsDb` annotation
#'   matching, unless `gffFile` is set.
#' @param gffFile `string` or `NULL`. *Advanced use only; not recommended.* By
#'   default, we recommend leaving this `NULL` for genomes that are supported on
#'   Ensembl. In this case, the row annotations ([rowRanges()]) will be obtained
#'   automatically from Ensembl by passing the `organism`, `genomeBuild`, and
#'   `ensemblRelease` arguments to AnnotationHub and ensembldb. For a genome
#'   that is not supported on Ensembl and/or AnnotationHub, a GFF/GTF (General
#'   Feature Format) file is required. Generally, we recommend using a GTF
#'   (GFFv2) file here over a GFF3 file if possible, although all GFF formats
#'   are supported. The function will internally generate transcript-to-gene
#'   mappings and construct a `GRanges` object containing the genomic ranges
#'   ([rowRanges()]).
#' @param ... Additional arguments, to be stashed in the [metadata()] slot.
#'
#' @return `bcbioSingleCell`.
#'
#' @seealso
#' - [SingleCellExperiment::SingleCellExperiment()].
#' - `.S4methods(class = "bcbioSingleCell")`.
#'
#' @examples
#' uploadDir <- system.file("extdata/indrops", package = "bcbioSingleCell")
#'
#' x <- bcbioSingleCell(uploadDir)
#' print(x)
#'
#' x <- bcbioSingleCell(
#'     uploadDir = uploadDir,
#'     sampleMetadataFile = file.path(uploadDir, "metadata.csv")
#' )
#' print(x)
bcbioSingleCell <- function(
    uploadDir,
    sampleMetadataFile = NULL,
    organism = NULL,
    ensemblRelease = NULL,
    genomeBuild = NULL,
    gffFile = NULL,
    transgeneNames = NULL,
    spikeNames = NULL,
    interestingGroups = "sampleName",
    ...
) {
    allSamples <- TRUE
    sampleData <- NULL

    # Legacy arguments ---------------------------------------------------------
    dots <- list(...)
    call <- match.call()
    # ensemblVersion
    if ("ensemblVersion" %in% names(call)) {
        warning("Use `ensemblRelease` instead of `ensemblVersion`.")
        ensemblRelease <- call[["ensemblVersion"]]
        dots[["ensemblVersion"]] <- NULL
    }
    # gtfFile
    if ("gtfFile" %in% names(call)) {
        warning("Use `gffFile` instead of `gtfFile`.")
        gffFile <- call[["gtfFile"]]
        dots[["gtfFile"]] <- NULL
    }
    # annotable
    if ("annotable" %in% names(call)) {
        stop("Use `gffFile` instead of `annotable`.")
    }
    rm(dots, call)

    # Assert checks ------------------------------------------------------------
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    assertIsAStringOrNULL(sampleMetadataFile)
    assertIsAStringOrNULL(organism)
    assertIsAnImplicitIntegerOrNULL(ensemblRelease)
    assertIsAStringOrNULL(genomeBuild)
    assertIsAStringOrNULL(gffFile)
    if (is_a_string(gffFile)) {
        assert_all_are_existing_files(gffFile)
    }
    assert_is_any_of(transgeneNames, c("character", "NULL"))
    assert_is_any_of(spikeNames, c("character", "NULL"))
    assert_is_character(interestingGroups)

    # Directory paths ----------------------------------------------------------
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)
    projectDir <- projectDir(uploadDir)
    sampleDirs <- sampleDirs(uploadDir)

    # Sequencing lanes ---------------------------------------------------------
    lanes <- detectLanes(sampleDirs)

    # Project summary YAML -----------------------------------------------------
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    yaml <- import(yamlFile)

    # bcbio run information ----------------------------------------------------
    dataVersions <-
        readDataVersions(file.path(projectDir, "data_versions.csv"))
    programVersions <-
        readProgramVersions(file.path(projectDir, "programs.txt"))
    log <-
        readLog(file.path(projectDir, "bcbio-nextgen.log"))
    commandsLog <-
        readLog(file.path(projectDir, "bcbio-nextgen-commands.log"))

    assert_is_tbl_df(dataVersions)
    assert_is_tbl_df(programVersions)
    assert_is_character(log)
    assert_is_character(commandsLog)

    cutoff <- getBarcodeCutoffFromCommandsLog(commandsLog)
    level <- getLevelFromCommandsLog(commandsLog)
    umiType <- getUMITypeFromCommandsLog(commandsLog)

    # Check to see if we're dealing with a multiplexed platform.
    multiplexed <- any(vapply(
        X = c("dropseq", "indrop"),
        FUN = function(pattern) {
            grepl(pattern = pattern, x = umiType)
        },
        FUN.VALUE = logical(1L)
    ))

    # User-defined sample metadata ---------------------------------------------
    if (is_a_string(sampleMetadataFile)) {
        sampleData <- readSampleData(sampleMetadataFile, lanes = lanes)

        # Allow sample selection by with this file.
        if (nrow(sampleData) < length(sampleDirs)) {
            message("Loading a subset of samples, defined by the metadata.")
            allSamples <- FALSE
            sampleDirs <- sampleDirs[rownames(sampleData)]
            message(paste(length(sampleDirs), "samples matched by metadata."))
        }

        # Error on incorrect reverse complement input.
        if ("sequence" %in% colnames(sampleData)) {
            sampleDirSequence <- str_extract(names(sampleDirs), "[ACGT]+$")
            if (identical(
                sort(sampleDirSequence),
                sort(as.character(sampleData[["sequence"]]))
            )) {
                stop(paste(
                    "It appears that the reverse complement sequence of the",
                    "i5 index barcodes were input into the sample metadata",
                    "`sequence` column. bcbio outputs the revcomp into the",
                    "sample directories, but the forward sequence should be",
                    "used in the R package."
                ))
            }
        }
    }

    # Unfiltered cellular barcode distributions --------------------------------
    cbList <- .import.bcbio.barcodes(sampleDirs)

    # Assays -------------------------------------------------------------------
    # Note that we're now allowing transcript-level counts (as of v0.99).
    counts <- .import.bcbio(sampleDirs)

    # Row data -----------------------------------------------------------------
    if (is_a_string(gffFile)) {
        rowRanges <- makeGRangesFromGFF(gffFile, level = level)
    } else if (is_a_string(organism)) {
        # Using AnnotationHub/ensembldb to obtain the annotations.
        message("Using `makeGRangesFromEnsembl()` for annotations.")
        rowRanges <- makeGRangesFromEnsembl(
            organism = organism,
            level = level,
            build = genomeBuild,
            release = ensemblRelease
        )
        if (is.null(genomeBuild)) {
            genomeBuild <- metadata(rowRanges)[["build"]]
        }
        if (is.null(ensemblRelease)) {
            ensemblRelease <- metadata(rowRanges)[["release"]]
        }
    } else {
        message("Unknown organism. Skipping annotations.")
        rowRanges <- emptyRanges(rownames(counts))
    }
    assert_is_all_of(rowRanges, "GRanges")
    assert_is_subset(rownames(counts), names(rowRanges))

    # Column data --------------------------------------------------------------
    # Automatic sample metadata.
    if (is.null(sampleData)) {
        if (isTRUE(multiplexed)) {
            # Multiplexed samples without user-defined metadata.
            message(paste0(
                "`sampleMetadataFile` is recommended for ",
                "multiplexed samples (e.g. ", umiType, ")."
            ))
            sampleData <- minimalSampleData(basename(sampleDirs))
        } else {
            sampleData <- readYAMLSampleData(yamlFile)
        }
    }
    assert_is_subset(rownames(sampleData), names(sampleDirs))

    # Always prefilter, removing very low quality cells with no UMIs or genes.
    colData <- metrics(
        object = counts,
        rowRanges = rowRanges,
        prefilter = TRUE
    )

    # Subset the counts to match the prefiltered metrics.
    assert_is_subset(rownames(colData), colnames(counts))
    counts <- counts[, rownames(colData), drop = FALSE]

    # Bind the `nCount` column into the colData. These are the number of counts
    # bcbio uses for initial filtering (minimum_barcode_depth in YAML).
    nCount <- .nCount(cbList)
    assert_is_subset(rownames(colData), names(nCount))
    colData[["nCount"]] <- nCount[rownames(colData)]

    # Join `sampleData` into cell-level `colData`.
    if (has_length(nrow(sampleData), n = 1L)) {
        colData[["sampleID"]] <- as.factor(rownames(sampleData))
    } else {
        colData[["sampleID"]] <- mapCellsToSamples(
            cells = rownames(colData),
            samples = rownames(sampleData)
        )
    }
    sampleData[["sampleID"]] <- as.factor(rownames(sampleData))
    colData <- left_join(
        x = colData,
        y = sampleData,
        by = "sampleID"
    )

    # Metadata -----------------------------------------------------------------
    # TODO Make this a function in bcbioBase.
    # Run date and template name.
    match <- str_match(
        string = basename(projectDir),
        pattern = projectDirPattern
    )
    runDate <- as.Date(match[[2L]])
    template <- match[[3L]]
    rm(match)

    # Interesting groups.
    interestingGroups <- camel(interestingGroups)
    assert_is_subset(interestingGroups, colnames(sampleData))

    metadata <- list(
        version = packageVersion,
        pipeline = "bcbio",
        level = level,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = as.character(sampleMetadataFile),
        interestingGroups = interestingGroups,
        organism = as.character(organism),
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease),
        umiType = umiType,
        allSamples = allSamples,
        lanes = lanes,
        # bcbio-specific -------------------------------------------------------
        projectDir = projectDir,
        template = template,
        runDate = runDate,
        yaml = yaml,
        gffFile = as.character(gffFile),
        tx2gene = tx2gene,
        dataVersions = dataVersions,
        programVersions = programVersions,
        bcbioLog = log,
        bcbioCommandsLog = commandsLog,
        cellularBarcodes = cbList,
        cellularBarcodeCutoff = cutoff,
        call = match.call()
    )

    # Return -------------------------------------------------------------------
    .new.bcbioSingleCell(
        assays = list(counts = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
}



.new.bcbioSingleCell <-  # nolint
    function(...) {
        new(
            Class = "bcbioSingleCell",
            makeSingleCellExperiment(...)
        )
    }



# CellRanger ===================================================================
# FIXME Can we parse the CellRanger `runDate` from the refData YAML?
# FIXME Allow this function to work if the user points at dir containing matrix.mtx.
# FIXME Add documentation about simple mode.

#' Read 10X Genomics Cell Ranger Single-Cell RNA-Seq Data
#'
#' Read [10x Genomics Cell Ranger](https://www.10xgenomics.com/software/)
#' counts into a `SingleCellExperiment` object.
#'
#' @section Directory structure:
#' Cell Ranger can vary in its output directory structure, but we're requiring a
#' single, consistent directory structure for all datasets, even those that only
#' contain a single sample:
#'
#' \preformatted{
#' file.path(
#'     "<uploadDir>",
#'     "<sampleName>",
#'     "outs",
#'     "filtered_gene_bc_matrices*",
#'     "<genomeBuild>",
#'     "matrix.mtx"
#' )
#' }
#'
#' @section Sample metadata:
#' A user-supplied sample metadata file defined by `sampleMetadataFile` is
#' required for multiplexed datasets. Otherwise this can be left `NULL`, and
#' minimal sample data will be used, based on the directory names.
#'
#' @section Reference data:
#' We strongly recommend supplying the corresponding reference data required for
#' Cell Ranger with the `refdataDir` argument. It will convert the gene
#' annotations defined in the GTF file into a `GRanges` object, which get
#' slotted in [rowRanges()]. Otherwise, the function will attempt to use the
#' most current annotations available from Ensembl, and some gene IDs may not
#' match, due to deprecation in the current Ensembl release.
#'
#' @family S4 Generators
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams bcbioSingleCell
#' @inheritParams general
#' @param uploadDir `string`. Path to Cell Ranger output directory. This
#'   directory path must contain `filtered_gene_bc_matrices*` as a child
#'   directory.
#' @param filtered `boolean`. Use filtered (recommended) or raw counts. Note
#'   that raw counts still contain only whitelisted cellular barcodes.
#' @param format `string`. Output format, either MatrixMarket ("`mtx`") or HDF5
#'   ("`hdf5`").
#' @param refdataDir `string` or `NULL`. Directory path to Cell Ranger reference
#'   annotation data.
#'
#' @return `CellRanger`.
#'
#' @examples
#' uploadDir <- system.file("extdata/cellranger", package = "bcbioSingleCell")
#' x <- CellRanger(uploadDir)
#' print(x)
CellRanger <- function(
    uploadDir,
    format = c("mtx", "hdf5"),
    filtered = TRUE,
    sampleMetadataFile = NULL,
    organism = NULL,
    ensemblRelease = NULL,
    genomeBuild = NULL,
    refdataDir = NULL,
    gffFile = NULL,
    transgeneNames = NULL,
    spikeNames = NULL,
    interestingGroups = "sampleName"
) {
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)
    format <- match.arg(format)
    assert_is_a_bool(filtered)
    assertIsAStringOrNULL(sampleMetadataFile)
    assertIsAStringOrNULL(organism)
    assertIsAnImplicitIntegerOrNULL(ensemblRelease)
    assertIsAStringOrNULL(genomeBuild)
    assertIsAStringOrNULL(refdataDir)
    if (is_a_string(refdataDir)) {
        assert_all_are_dirs(refdataDir)
        refdataDir <- normalizePath(refdataDir, winslash = "/", mustWork = TRUE)
    }
    assertIsAStringOrNULL(gffFile)
    assert_is_any_of(transgeneNames, c("character", "NULL"))
    assert_is_any_of(spikeNames, c("character", "NULL"))
    assert_is_character(interestingGroups)

    pipeline <- "cellranger"
    level <- "genes"
    umiType <- "chromium"

    # Sample files -------------------------------------------------------------
    sampleFiles <- .sampleFiles.cellranger(uploadDir)

    # Sequencing lanes ---------------------------------------------------------
    lanes <- detectLanes(sampleFiles)

    # Sample metadata ----------------------------------------------------------
    allSamples <- TRUE
    sampleData <- NULL

    if (is_a_string(sampleMetadataFile)) {
        sampleData <- readSampleData(sampleMetadataFile)
        # Allow sample selection by with this file.
        if (nrow(sampleData) < length(sampleFiles)) {
            message("Loading a subset of samples, defined by the metadata.")
            allSamples <- FALSE
            sampleFiles <- sampleFiles[rownames(sampleData)]
            message(paste(length(sampleFiles), "samples matched by metadata."))
        }
    }

    # Counts -------------------------------------------------------------------
    counts <- .import.cellranger(sampleFiles)

    # Row data -----------------------------------------------------------------
    refJSON <- NULL

    # Prepare gene annotations as GRanges.
    if (is_a_string(refdataDir)) {
        message("Using 10X Genomics reference data for annotations.")
        message(paste("refdataDir:", refdataDir))
        # JSON data.
        refJSONFile <- file.path(refdataDir, "reference.json")
        assert_all_are_existing_files(refJSONFile)
        refJSON <- import(refJSONFile)
        # Get the genome build from JSON metadata.
        genomeBuild <- unlist(refJSON[["genomes"]])
        assert_is_a_string(genomeBuild)
        # Convert the GTF file to GRanges.
        gffFile <- file.path(refdataDir, "genes", "genes.gtf")
        assert_is_a_string(gffFile)
        rowRanges <- makeGRangesFromGFF(gffFile)
        # Get the Ensembl version from the GTF file name.
        # Example: "Homo_sapiens.GRCh37.82.filtered.gtf"
        ensemblRelease <- gffFile %>%
            str_split("\\.", simplify = TRUE) %>%
            .[1L, 3L] %>%
            as.integer()
    } else if (is_a_string(gffFile)) {
        rowRanges <- makeGRangesFromGFF(gffFile, level = "genes")
    } else if (is_a_string(organism)) {
        # CellRanger uses Ensembl refdata internally. Here we're fetching the
        # annotations with AnnotationHub rather than pulling from the GTF file
        # in the refdata directory. It will also drop genes that are now dead in
        # the current Ensembl release. Don't warn about old Ensembl release
        # version.
        message("Using `makeGRangesFromEnsembl()` for annotations.")
        rowRanges <- makeGRangesFromEnsembl(
            organism = organism,
            level = level,
            build = genomeBuild,
            release = ensemblRelease
        )
        if (is.null(genomeBuild)) {
            genomeBuild <- metadata(rowRanges)[["build"]]
        }
        if (is.null(ensemblRelease)) {
            ensemblRelease <- metadata(rowRanges)[["release"]]
        }
    } else {
        message("Unknown organism. Skipping annotations.")
        rowRanges <- emptyRanges(rownames(counts))
    }
    assert_is_all_of(rowRanges, "GRanges")

    # Column data --------------------------------------------------------------
    # Automatic sample metadata.
    if (is.null(sampleData)) {
        # Define the grep pattern to use for sample ID extraction.
        pattern <- "^(.+)_[ACGT]+$"
        if (all(grepl(pattern, colnames(counts)))) {
            match <- str_match(
                string = colnames(counts),
                pattern = pattern
            )
            samples <- unique(match[, 2L, drop = TRUE])
        } else if (has_length(sampleFiles, n = 1L)) {
            samples <- names(sampleFiles)
        }
        sampleData <- minimalSampleData(samples)
    }

    # Always prefilter, removing very low quality cells with no UMIs or genes.
    colData <- metrics(
        object = counts,
        rowRanges = rowRanges,
        prefilter = TRUE
    )

    # Subset the counts to match the prefiltered metrics.
    assert_is_subset(rownames(colData), colnames(counts))
    counts <- counts[, rownames(colData), drop = FALSE]

    # Join `sampleData` into cell-level `colData`.
    if (has_length(nrow(sampleData), n = 1L)) {
        colData[["sampleID"]] <- as.factor(rownames(sampleData))
    } else {
        colData[["sampleID"]] <- mapCellsToSamples(
            cells = rownames(colData),
            samples = rownames(sampleData)
        )
    }

    # Metadata -----------------------------------------------------------------
    # Interesting groups.
    interestingGroups <- camel(interestingGroups)
    assert_is_subset(interestingGroups, colnames(sampleData))

    metadata <- list(
        version = packageVersion,
        pipeline = pipeline,
        level = level,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = as.character(sampleMetadataFile),
        interestingGroups = interestingGroups,
        organism = organism,
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease),
        umiType = umiType,
        allSamples = allSamples,
        lanes = lanes,
        # cellranger-specific --------------------------------------------------
        refdataDir = refdataDir,
        refJSON = refJSON,
        call = match.call()
    )

    # Return -------------------------------------------------------------------
    .new.CellRanger(
        assays = list(counts = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
}



.new.CellRanger <-  # nolint
    function(...) {
        new(
            Class = "CellRanger",
            makeSingleCellExperiment(...)
        )
    }
