#' Read bcbio Single-Cell RNA-Seq Data
#'
#' @section Remote Data:
#'
#' When working in RStudio, we recommend connecting to the bcbio-nextgen run
#' directory as a remote connection over
#' [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
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
    uploadDir <- realpath(uploadDir)
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
    assert_is_all_of(dataVersions, "DataFrame")

    programVersions <-
        readProgramVersions(file.path(projectDir, "programs.txt"))
    assert_is_all_of(programVersions, "DataFrame")

    log <-
        import(file.path(projectDir, "bcbio-nextgen.log"))
    assert_is_character(log)

    commands <-
        import(file.path(projectDir, "bcbio-nextgen-commands.log"))
    assert_is_character(commands)

    cutoff <- getBarcodeCutoffFromCommands(commands)
    level <- getLevelFromCommands(commands)
    umiType <- getUMITypeFromCommands(commands)

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
            genomeBuild = genomeBuild,
            release = ensemblRelease
        )
        if (is.null(genomeBuild)) {
            genomeBuild <- metadata(rowRanges)[["genomeBuild"]]
        }
        if (is.null(ensemblRelease)) {
            ensemblRelease <- metadata(rowRanges)[["ensemblRelease"]]
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
            sampleData <- getSampleDataFromYAML(yaml)
        }
    }
    assert_is_subset(rownames(sampleData), names(sampleDirs))

    # Always prefilter, removing very low quality cells with no UMIs or genes.
    colData <- metrics.matrix(
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
    assert_is_integer(nCount)
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
        x = as_tibble(colData, rownames = "rowname"),
        y = as_tibble(sampleData, rownames = NULL),
        by = "sampleID"
    )
    colData <- as(colData, "DataFrame")

    # Metadata -----------------------------------------------------------------
    runDate <- runDate(projectDir)

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
        runDate = runDate,
        yaml = yaml,
        gffFile = as.character(gffFile),
        dataVersions = dataVersions,
        programVersions = programVersions,
        bcbioLog = log,
        bcbioCommandsLog = commands,
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
