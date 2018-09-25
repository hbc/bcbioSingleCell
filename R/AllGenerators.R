# bcbioSingleCell ==============================================================
#' Read bcbio Single-Cell RNA-Seq Data
#'
#' @section Remote Data:
#'
#' When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @note `bcbioSingleCell` extended `SummarizedExperiment` prior to v0.1.0,
#'   where we migrated to `SingleCellExperiment`.
#'
#' @family S4 Generators
#' @author Michael Steinbaugh, Rory Kirchner
#' @export
#'
#' @inheritParams basejump::makeSummarizedExperiment
#' @inheritParams general
#' @param uploadDir `string`. Path to final upload directory. This path is set
#'   when running "`bcbio_nextgen -w template`".
#' @param organism `string` or `NULL`. Organism name. Use the full Latin name
#'   (e.g. "Homo sapiens"), since this will be input downstream to AnnotationHub
#'   and ensembldb, unless `gffFile` is set. If left `NULL` (*not recommended*),
#'   the function call will skip loading gene-level annotations into
#'   [rowRanges()]. This can be useful for poorly annotation genomes or
#'   experiments involving multiple genomes.
#' @param sampleMetadataFile `string` or `NULL`. Sample barcode metadata file.
#'   Optional for runs with demultiplixed index barcodes (e.g. SureCell) and can
#'   be left `NULL`, but otherwise required for runs with multipliexed FASTQs
#'   containing multiple index barcodes (e.g. inDrop).
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
#' x <- bcbioSingleCell(
#'     uploadDir = uploadDir,
#'     organism = "Homo sapiens",
#'     sampleMetadataFile = file.path(uploadDir, "metadata.csv"),
#'     ensemblRelease = 87L
#' )
#' print(x)
bcbioSingleCell <- function(
    uploadDir,
    organism = NULL,
    sampleMetadataFile = NULL,
    interestingGroups = "sampleName",
    ensemblRelease = NULL,
    genomeBuild = NULL,
    transgeneNames = NULL,
    spikeNames = NULL,
    gffFile = NULL,
    ...
) {
    dots <- list(...)
    pipeline <- "bcbio"

    # Legacy arguments ---------------------------------------------------------
    call <- match.call(expand.dots = TRUE)
    # annotable
    if ("annotable" %in% names(call)) {
        stop("Use `gffFile` instead of `annotable`")
    }
    # ensemblVersion
    if ("ensemblVersion" %in% names(call)) {
        warning("Use `ensemblRelease` instead of `ensemblVersion`")
        ensemblRelease <- call[["ensemblVersion"]]
        dots[["ensemblVersion"]] <- NULL
    }
    # gtfFile
    if ("gtfFile" %in% names(call)) {
        warning("Use `gffFile` instead of `gtfFile`")
        gffFile <- call[["gtfFile"]]
        dots[["gtfFile"]] <- NULL
    }
    # organism
    if (!"organism" %in% names(call)) {
        message("`organism` is now recommended, to acquire gene annotations")
    }
    dots <- Filter(Negate(is.null), dots)

    # Assert checks ------------------------------------------------------------
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    assertIsAStringOrNULL(sampleMetadataFile)
    assert_is_character(interestingGroups)
    assertIsAStringOrNULL(organism)
    assertIsAnImplicitIntegerOrNULL(ensemblRelease)
    assertIsAStringOrNULL(genomeBuild)
    assert_is_any_of(transgeneNames, c("character", "NULL"))
    assert_is_any_of(spikeNames, c("character", "NULL"))
    assertIsAStringOrNULL(gffFile)
    if (is_a_string(gffFile)) {
        assert_all_are_existing_files(gffFile)
    }

    # Directory paths ----------------------------------------------------------
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)
    projectDir <- projectDir(uploadDir)
    sampleDirs <- sampleDirs(uploadDir)

    # Run date and template name -----------------------------------------------
    # Get run date and template name from project directory.
    # This information will be stashed in `metadata()`.
    match <- str_match(
        string = basename(projectDir),
        pattern = projectDirPattern
    )
    runDate <- as.Date(match[[2L]])
    template <- match[[3L]]
    rm(match)

    # Sequencing lanes ---------------------------------------------------------
    if (any(grepl(lanePattern, sampleDirs))) {
        lanes <- str_match(names(sampleDirs), lanePattern) %>%
            .[, 2L] %>%
            unique() %>%
            length()
        message(paste(
            lanes, "sequencing lane detected", "(technical replicates)"
        ))
    } else {
        lanes <- 1L
    }

    # Project summary YAML -----------------------------------------------------
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    yaml <- readYAML(yamlFile)

    # bcbio run information ----------------------------------------------------
    dataVersions <- readDataVersions(
        file = file.path(projectDir, "data_versions.csv")
    )
    assert_is_tbl_df(dataVersions)

    programVersions <- readProgramVersions(
        file = file.path(projectDir, "programs.txt")
    )
    assert_is_tbl_df(programVersions)

    bcbioLog <- readLog(
        file = file.path(projectDir, "bcbio-nextgen.log")
    )
    assert_is_character(bcbioLog)

    bcbioCommandsLog <- readLog(
        file = file.path(projectDir, "bcbio-nextgen-commands.log")
    )
    assert_is_character(bcbioCommandsLog)

    # Cellular barcode cutoff --------------------------------------------------
    cellularBarcodeCutoffPattern <- "--cb_cutoff (\\d+)"
    assert_any_are_matching_regex(
        x = bcbioCommandsLog,
        pattern = cellularBarcodeCutoffPattern
    )
    match <- str_match(
        string = bcbioCommandsLog,
        pattern = cellularBarcodeCutoffPattern
    )
    cellularBarcodeCutoff <- match %>%
        .[, 2L] %>%
        na.omit() %>%
        unique() %>%
        as.integer()
    assert_is_an_integer(cellularBarcodeCutoff)

    message(paste(
        cellularBarcodeCutoff,
        "reads per cellular barcode cutoff detected"
    ))

    # Detect gene or transcript-level output -----------------------------------
    genemapPattern <- "--genemap (.+)-tx2gene.tsv"
    if (any(grepl(genemapPattern, bcbioCommandsLog))) {
        level <- "genes"
    } else {
        level <- "transcripts"
    }

    # Molecular barcode (UMI) type ---------------------------------------------
    umiPattern <- "fastqtransform.*/(.*)\\.json"
    assert_any_are_matching_regex(bcbioCommandsLog, umiPattern)
    umiType <- str_match(bcbioCommandsLog, umiPattern) %>%
        .[, 2L] %>%
        na.omit() %>%
        unique() %>%
        gsub(
            pattern = "-transform",
            replacement = "",
            x = .
        )
    assert_is_a_string(umiType)
    message(paste("UMI type:", umiType))

    # Sample metadata ----------------------------------------------------------
    # External file required for inDrop
    if (grepl("indrop", umiType) && is.null(sampleMetadataFile)) {
        stop(paste(
            "inDrop samples require `sampleMetadataFile`",
            "containing the index barcode sequences"
        ))
    }

    if (is_a_string(sampleMetadataFile)) {
        sampleData <- readSampleData(sampleMetadataFile)
    } else {
        sampleData <- readYAMLSampleData(yamlFile)
    }

    # Check for incorrect reverse complement input
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

    assert_is_subset(rownames(sampleData), names(sampleDirs))
    sampleData <- sanitizeSampleData(sampleData)

    # Interesting groups -------------------------------------------------------
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups)
    assert_is_subset(interestingGroups, colnames(sampleData))

    # Subset sample directories by metadata ------------------------------------
    # Check to see if a subset of samples is requested via the metadata file.
    # This matches by the reverse complement sequence of the index barcode.
    if (nrow(sampleData) < length(sampleDirs)) {
        message("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs[rownames(sampleData)]
        message(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Assays -------------------------------------------------------------------
    message(paste("Reading counts as", level))
    counts <- .readCounts(
        sampleDirs = sampleDirs,
        pipeline = pipeline,
        format = "mtx"
    )

    # Require transcript to gene conversion (legacy) ---------------------------
    if (level == "transcripts") {
        message("Converting transcripts to genes")

        # GFF file is required
        if (!is_a_string(gffFile)) {
            stop("GFF is required to convert transcripts to genes")
        }

        # Note that this data.frame won't contain transcript versions
        tx2gene <- makeTx2geneFromGFF(gffFile)

        # Add spike-ins to tx2gene, if necessary
        if (is.character(isSpike)) {
            assert_are_disjoint_sets(rownames(tx2gene), isSpike)
            spike <- data.frame(
                transcriptID = isSpike,
                geneID = isSpike,
                row.names = isSpike,
                stringsAsFactors = FALSE
            )
            tx2gene <- rbind(spike, tx2gene)
        }

        # Ensure Ensembl transcript versions are removed
        counts <- stripTranscriptVersions(counts)

        # Resize the tx2gene to match the matrix rownames
        assert_is_subset(rownames(counts), rownames(tx2gene))
        tx2gene <- tx2gene[rownames(counts), , drop = FALSE]

        # Now we're ready to assign `geneID` and aggregate
        rownames(counts) <- tx2gene[["geneID"]]
        counts <- aggregate.Matrix(
            x = counts,
            groupings = rownames(counts),
            fun = "sum"
        )
        level <- "genes"
    } else {
        tx2gene <- NULL
    }

    # Unfiltered cellular barcode distributions --------------------------------
    cbList <- .cellularBarcodesList(sampleDirs)
    cbData <- .bindCellularBarcodes(cbList)

    # Row data -----------------------------------------------------------------
    if (is_a_string(gffFile)) {
        rowRanges <- makeGRangesFromGFF(gffFile, level = "genes")
    } else if (is_a_string(organism)) {
        # Using AnnotationHub/ensembldb to obtain the annotations.
        message("Using `makeGRangesFromEnsembl()` for annotations")
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
        rowRanges <- emptyRanges(rownames(counts))
    }
    assert_is_all_of(rowRanges, "GRanges")
    assert_is_subset(rownames(counts), names(rowRanges))

    # Column data --------------------------------------------------------------
    # Always prefilter, removing very low quality cells with no UMIs or genes
    metrics <- .calculateMetrics(
        object = counts,
        rowRanges = rowRanges,
        prefilter = TRUE
    )

    # Subset the counts to match the prefiltered metrics
    counts <- counts[, rownames(metrics), drop = FALSE]

    colData <- as(metrics, "DataFrame")
    colData[["cellID"]] <- rownames(colData)
    cell2sample <- mapCellsToSamples(
        cells = rownames(colData),
        samples = rownames(sampleData)
    )
    colData[["sampleID"]] <- cell2sample
    sampleData[["sampleID"]] <- rownames(sampleData)
    colData <- merge(
        x = colData,
        y = sampleData,
        by = "sampleID",
        all.x = TRUE
    )
    rownames(colData) <- colData[["cellID"]]
    colData[["cellID"]] <- NULL
    sampleData[["sampleID"]] <- NULL

    # Bind the `nCount` column into the colData
    nCount <- cbData[rownames(colData), "nCount", drop = FALSE]
    colData <- cbind(nCount, colData)

    # Metadata -----------------------------------------------------------------
    metadata <- list(
        version = packageVersion,
        pipeline = pipeline,
        level = level,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = as.character(sampleMetadataFile),
        interestingGroups = interestingGroups,
        organism = as.character(organism),
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease),
        cell2sample = as.factor(cell2sample),
        umiType = umiType,
        allSamples = allSamples,
        # bcbio pipeline-specific ----------------------------------------------
        projectDir = projectDir,
        template = template,
        runDate = runDate,
        yaml = yaml,
        gffFile = as.character(gffFile),
        tx2gene = tx2gene,
        dataVersions = dataVersions,
        programVersions = programVersions,
        bcbioLog = bcbioLog,
        bcbioCommandsLog = bcbioCommandsLog,
        cellularBarcodes = cbList,
        cellularBarcodeCutoff = cellularBarcodeCutoff,
        call = match.call()
    )
    # Add user-defined custom metadata, if specified
    if (length(dots)) {
        assert_are_disjoint_sets(names(metadata), names(dots))
        metadata <- c(metadata, dots)
    }

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



# CellRanger ===================================================================
# FIXME Can we parse the CellRanger `runDate` from the YAML?

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
    organism = NULL,
    sampleMetadataFile = NULL,
    interestingGroups = "sampleName",
    refdataDir = NULL,
    gffFile = NULL,
    transgeneNames = NULL,
    spikeNames = NULL,
    ...
) {
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)
    format <- match.arg(format)
    assert_is_a_bool(filtered)
    assertIsAStringOrNULL(organism)
    assertIsAStringOrNULL(sampleMetadataFile)
    assertIsAStringOrNULL(refdataDir)
    if (is_a_string(refdataDir)) {
        assert_all_are_dirs(refdataDir)
        refdataDir <- normalizePath(refdataDir, winslash = "/", mustWork = TRUE)
    }
    assertIsAStringOrNULL(gffFile)
    assert_is_character(interestingGroups)
    assert_is_any_of(transgeneNames, c("character", "NULL"))
    assert_is_any_of(spikeNames, c("character", "NULL"))
    dots <- list(...)
    pipeline <- "cellranger"
    level <- "genes"
    umiType <- "chromium"

    # Legacy arguments ---------------------------------------------------------
    # annotable
    stopifnot(!"annotable" %in% names(call))
    dots <- Filter(Negate(is.null), dots)

    # Sample directories -------------------------------------------------------
    dirs <- list.dirs(uploadDir, recursive = FALSE)
    assert_is_non_empty(dirs)
    # Sample subdirectories must contain `outs/` directory
    hasOuts <- vapply(
        X = dirs,
        FUN = function(dir) {
            dir.exists(file.path(dir, "outs"))
        },
        FUN.VALUE = logical(1L)
    )
    sampleDirs <- dirs[hasOuts]
    assert_is_non_empty(sampleDirs)
    names(sampleDirs) <- makeNames(basename(sampleDirs), unique = TRUE)
    message(paste(length(sampleDirs), "sample(s) detected"))

    # Sample metadata ----------------------------------------------------------
    if (is_a_string(sampleMetadataFile)) {
        sampleData <- readSampleData(sampleMetadataFile)
    } else {
        sampleData <- minimalSampleData(basename(sampleDirs))
    }

    # Interesting groups -------------------------------------------------------
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups)
    assert_is_subset(interestingGroups, colnames(sampleData))

    # Subset sample directories by metadata ------------------------------------
    if (nrow(sampleData) < length(sampleDirs)) {
        message("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs[rownames(sampleData)]
        message(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Counts -------------------------------------------------------------------
    # This step can be slow over sshfs, recommend running on an HPC
    message("Reading counts at gene level")
    counts <- .readCounts(
        sampleDirs = sampleDirs,
        pipeline = pipeline,
        format = format,
        filtered = filtered
    )

    # Multiplexed sample check -------------------------------------------------
    # Check to see if multiplexed samples are present and require metadata
    multiplexedPattern <- "^(.+)_(\\d+)_([ACGT]+)$"
    if (any(grepl(multiplexedPattern, colnames(counts)))) {
        message("Multiplexed (aggregated) samples detected")
        # Prepare data.frame of barcode mappings
        # Example:
        # cellID: cellranger_AAACCTGGTTTACTCT_1
        # description: cellranger
        # barcode: AAACCTGGTTTACTCT
        # index: 1
        cellMap <- str_match(
            string = colnames(counts),
            pattern = "^(.+)_(\\d+)_([ACGT]+)$"
        ) %>%
            as.data.frame() %>%
            set_colnames(c(
                "cellID",
                "description",
                "index",
                "barcode"
            )) %>%
            mutate_all(as.factor)

        # Check for single sample and fix sampleData automatically if necessary
        if (
            identical(levels(cellMap[["index"]]), "1") &&
            !"index" %in% colnames(sampleData)
        ) {
            sampleData[["index"]] <- factor("1")
            rownames(sampleData) <- paste0(rownames(sampleData), "_1")
        }

        # Require user defined metadata
        if (!"index" %in% colnames(sampleData)) {
            stop(paste(
                "`index` column must be defined using",
                "`sampleMetadataFile` for multiplexed samples"
            ))
        }
    }

    # Row data -----------------------------------------------------------------
    refJSON <- NULL
    genomeBuild <- NULL
    ensemblRelease <- NULL

    # Prepare gene annotations as GRanges
    if (is_a_string(refdataDir)) {
        message("Using 10X Genomics reference data for gene annotations")
        message(paste("refdataDir:", refdataDir))
        # JSON data
        refJSONFile <- file.path(refdataDir, "reference.json")
        assert_all_are_existing_files(refJSONFile)
        refJSON <- readJSON(refJSONFile)
        # Get the genome build from JSON metadata
        genomeBuild <- unlist(refJSON[["genomes"]])
        assert_is_a_string(genomeBuild)
        # Convert the GTF file to GRanges
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
        # Using AnnotationHub/ensembldb to obtain the annotations.
        message("Using `makeGRangesFromEnsembl()` for annotations")
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
        rowRanges <- emptyRanges(rownames(counts))
    }
    assert_is_all_of(rowRanges, "GRanges")

    # Column data --------------------------------------------------------------
    # Always prefilter, removing very low quality cells with no UMIs or genes
    metrics <- .calculateMetrics(
        object = counts,
        rowRanges = rowRanges,
        prefilter = TRUE
    )

    # Subset the counts to match the prefiltered metrics
    counts <- counts[, rownames(metrics), drop = FALSE]

    colData <- as(metrics, "DataFrame")
    colData[["cellID"]] <- rownames(colData)
    cell2sample <- mapCellsToSamples(
        cells = rownames(colData),
        samples = rownames(sampleData)
    )
    colData[["sampleID"]] <- cell2sample
    sampleData[["sampleID"]] <- rownames(sampleData)
    colData <- merge(
        x = colData,
        y = sampleData,
        by = "sampleID",
        all.x = TRUE
    )
    rownames(colData) <- colData[["cellID"]]
    colData[["cellID"]] <- NULL
    sampleData[["sampleID"]] <- NULL

    # Metadata -----------------------------------------------------------------
    metadata <- list(
        version = packageVersion,
        pipeline = pipeline,
        level = level,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = as.character(sampleMetadataFile),
        interestingGroups = interestingGroups,
        cell2sample = as.factor(cell2sample),
        organism = organism,
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease),
        umiType = umiType,
        allSamples = allSamples,
        # cellranger pipeline-specific -----------------------------------------
        refdataDir = refdataDir,
        refJSON = refJSON,
        call = match.call()
    )
    # Add user-defined custom metadata, if specified
    if (length(dots)) {
        assert_are_disjoint_sets(names(metadata), names(dots))
        metadata <- c(metadata, dots)
    }

    # Return -------------------------------------------------------------------
    sce <- .new.SingleCellExperiment(
        assays = list(counts = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
    new("CellRanger", sce)
}
