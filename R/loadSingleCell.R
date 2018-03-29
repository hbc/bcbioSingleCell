#' Load bcbio Single-Cell RNA-Seq Data
#'
#' @note When working in RStudio, we recommend connecting to the bcbio-nextgen
#'   run directory as a remote connection over
#'   [sshfs](https://github.com/osxfuse/osxfuse/wiki/SSHFS).
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#' @param uploadDir Path to final upload directory. This path is set when
#'   running `bcbio_nextgen -w template`.
#' @param sampleMetadataFile Sample barcode metadata file. Optional for runs
#'   with demultiplixed index barcodes (e.g. SureCell), but otherwise required
#'   for runs with multipliexed FASTQs containing multiple index barcodes (e.g.
#'   inDrop). Consult the GitHub repo for examples and additional information.
#' @param organism Organism name. Use the full latin name (e.g.
#'   "Homo sapiens"), since this will be input downstream to
#'   AnnotationHub and ensembldb, unless `gffFile` is set.
#' @param genomeBuild *Optional.* Ensembl genome build name (e.g. "GRCh38").
#'   This will be passed to AnnotationHub for `EnsDb` annotation matching,
#'   unless `gffFile` is set.
#' @param ensemblRelease *Optional.* Ensembl release version. If `NULL`,
#'   defaults to current release, and does not typically need to be
#'   user-defined. Passed to AnnotationHub for `EnsDb` annotation matching,
#'   unless `gffFile` is set.
#' @param isSpike *Optional.* Gene names corresponding to FASTA spike-in
#'   sequences (e.g. ERCCs, EGFP, TDTOMATO).
#' @param gffFile *Optional, not recommended.* By default, we recommend leaving
#'   this `NULL` for genomes that are supported on Ensembl. In this case, the
#'   row annotations ([rowRanges()]) will be obtained automatically from Ensembl
#'   by passing the `organism`, `genomeBuild`, and `ensemblRelease` arguments to
#'   AnnotationHub and ensembldb. For a genome that is not supported on Ensembl
#'   and/or AnnotationHub, a GFF/GTF (General Feature Format) file is required.
#'   Generally, we recommend using a GTF (GFFv2) file here over a GFF3 file if
#'   possible, although all GFF formats are supported. The function will
#'   internally generate a `TxDb` containing transcript-to-gene mappings and
#'   construct a `GRanges` object containing the genomic ranges ([rowRanges()]).
#' @param prefilter Prefilter counts prior to quality control analysis.
#' @param ... Additional arguments, to be stashed in the [metadata()] slot.
#'
#' @return `bcbioSingleCell`.
#' @export
#'
#' @examples
#' uploadDir <- system.file("extdata/indrop", package = "bcbioSingleCell")
#' loadSingleCell(
#'     uploadDir = uploadDir,
#'     sampleMetadataFile = file.path(uploadDir, "metadata.csv"),
#'     organism = "Homo sapiens"
#' )
loadSingleCell <- function(
    uploadDir,
    sampleMetadataFile,
    interestingGroups = "sampleName",
    organism,
    genomeBuild = NULL,
    ensemblRelease = NULL,
    isSpike = NULL,
    gffFile = NULL,
    prefilter = TRUE,
    ...
) {
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)
    if (!missing(sampleMetadataFile)) {
        assert_is_a_string(sampleMetadataFile)
        # Allow for metadata from URL, so don't check if exists here
    }
    assert_is_character(interestingGroups)
    assert_is_a_string(organism)
    assertIsAStringOrNULL(genomeBuild)
    assertIsAnImplicitIntegerOrNULL(ensemblRelease)
    assertIsCharacterOrNULL(isSpike)
    assertIsAStringOrNULL(gffFile)
    if (is_a_string(gffFile)) {
        assert_all_are_existing_files(gffFile)
    }
    assert_is_a_bool(prefilter)
    dots <- list(...)
    pipeline <- "bcbio"

    # Legacy arguments =========================================================
    call <- match.call(expand.dots = TRUE)
    # annotable
    if ("annotable" %in% names(call)) {
        abort("Use `gffFile` instead of `annotable`")
    }
    # ensemblVersion
    if ("ensemblVersion" %in% names(call)) {
        warn("Use `ensemblRelease` instead of `ensemblVersion`")
        ensemblRelease <- call[["ensemblVersion"]]
        dots[["ensemblVersion"]] <- NULL
    }
    # gtfFile
    if ("gtfFile" %in% names(call)) {
        warn("Use `gffFile` instead of `gtfFile`")
        gffFile <- call[["gtfFile"]]
        dots[["gtfFile"]] <- NULL
    }
    dots <- Filter(Negate(is.null), dots)

    # Directory paths ==========================================================
    projectDir <- dir(
        uploadDir,
        pattern = bcbioBase::projectDirPattern,
        full.names = FALSE,
        recursive = FALSE
    )
    assert_is_a_string(projectDir)
    inform(projectDir)
    match <- str_match(projectDir, bcbioBase::projectDirPattern)
    runDate <- as.Date(match[[2L]])
    template <- match[[3L]]
    projectDir <- file.path(uploadDir, projectDir)
    sampleDirs <- .sampleDirs(uploadDir, pipeline = pipeline)

    # Sequencing lanes =========================================================
    if (any(grepl(bcbioBase::lanePattern, sampleDirs))) {
        lanes <- str_match(names(sampleDirs), bcbioBase::lanePattern) %>%
            .[, 2L] %>%
            unique() %>%
            length()
        inform(paste(
            lanes, "sequencing lane detected", "(technical replicates)"
        ))
    } else {
        lanes <- 1L
    }

    # Project summary YAML =====================================================
    yamlFile <- file.path(projectDir, "project-summary.yaml")
    yaml <- readYAML(yamlFile)

    # bcbio run information ====================================================
    dataVersions <- readDataVersions(
        file = file.path(projectDir, "data_versions.csv")
    )
    assert_is_tbl_df(dataVersions)

    programVersions <- readProgramVersions(
        file = file.path(projectDir, "programs.txt")
    )
    assert_is_tbl_df(programVersions)

    bcbioLog <- readLogFile(
        file = file.path(projectDir, "bcbio-nextgen.log")
    )
    assert_is_character(bcbioLog)

    bcbioCommandsLog <- readLogFile(
        file = file.path(projectDir, "bcbio-nextgen-commands.log")
    )
    assert_is_character(bcbioCommandsLog)

    # Cellular barcode cutoff ==================================================
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

    inform(paste(
        cellularBarcodeCutoff,
        "reads per cellular barcode cutoff detected"
    ))

    # Detect gene or transcript-level output ===================================
    genemapPattern <- "--genemap (.+)-tx2gene.tsv"
    if (any(grepl(genemapPattern, bcbioCommandsLog))) {
        level <- "genes"
    } else {
        level <- "transcripts"
    }

    # Molecular barcode (UMI) type =============================================
    umiPattern <- "/umis/([a-z0-9\\-]+)\\.json"
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
    inform(paste("UMI type:", umiType))

    # Sample metadata ==========================================================
    if (is_a_string(sampleMetadataFile)) {
        sampleData <- readSampleMetadataFile(sampleMetadataFile)
    } else {
        sampleData <- sampleYAMLMetadata(yaml)
    }
    # Check for reverse complement input
    if ("sequence" %in% colnames(sampleData)) {
        sampleDirSequence <- str_extract(names(sampleDirs), "[ACGT]+$")
        if (identical(
            sort(sampleDirSequence),
            sort(as.character(sampleData[["sequence"]]))
        )) {
            abort(paste(
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

    # Interesting groups =======================================================
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    assertFormalInterestingGroups(sampleData, interestingGroups)

    # Subset sample directories by metadata ====================================
    # Check to see if a subset of samples is requested via the metadata file.
    # This matches by the reverse complement sequence of the index barcode.
    if (nrow(sampleData) < length(sampleDirs)) {
        inform("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs[rownames(sampleData)]
        inform(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Row data =================================================================
    tx2gene <- NULL
    if (is_a_string(gffFile)) {
        rowRanges <- rowRangesFromGFF(gffFile, level = level)
        rowRangesMetadata <- NULL
        if (level == "transcripts") {
            tx2gene <- tx2geneFromGFF(gffFile)
        }
    } else {
        # ah = AnnotationHub
        ah <- ensembl(
            organism = organism,
            format = "genes",
            genomeBuild = genomeBuild,
            release = ensemblRelease,
            return = "GRanges",
            metadata = TRUE
        )
        assert_is_list(ah)
        assert_are_identical(names(ah), c("data", "metadata"))
        rowRanges <- ah[["data"]]
        assert_is_all_of(rowRanges, "GRanges")
        rowRangesMetadata <- ah[["metadata"]]
        assert_is_data.frame(rowRangesMetadata)
        if (level == "transcripts") {
            tx2gene <- ensembl(
                organism = organism,
                format = "tx2gene",
                genomeBuild = genomeBuild,
                release = ensemblRelease
            )
        }
    }

    # Require gene-to-symbol mappings
    assert_is_subset(
        x = c("geneID", "geneName"),
        y = names(mcols(rowRanges))
    )

    rowData <- as.data.frame(rowRanges)
    rownames(rowData) <- names(rowRanges)

    # Assays ===================================================================
    inform(paste("Reading counts as", level))
    countsList <- .sparseCountsList(
        sampleDirs = sampleDirs,
        pipeline = pipeline,
        umiType = umiType
    )
    # Ensure samples with empty matrices (`NULL`) are filtered
    countsList <- Filter(Negate(is.null), countsList)
    counts <- do.call(cBind, countsList)

    # Transcript to gene level counts (legacy)
    # Now recommended to provide GTF file during the bcbio run instead
    if (level == "transcript") {
        inform("Converting transcripts to genes")
        assert_is_subset(rownames(counts), rownames(tx2gene))
        # Resize the tx2gene to match the matrix rownames
        tx2gene <- tx2gene[rownames(counts), , drop = FALSE]
        assert_are_identical(rownames(counts), rownames(tx2gene))
        rownames(counts) <- tx2gene[["geneID"]]
        counts <- aggregate.Matrix(
            x = counts,
            groupings = rownames(counts),
            fun = "sum"
        )
    }

    # Unfiltered cellular barcode distributions ================================
    cbList <- .cellularBarcodesList(sampleDirs)
    cbData <- .bindCellularBarcodes(cbList)

    # Column data ==============================================================
    colData <- calculateMetrics(
        object = counts,
        rowData = rowData,
        prefilter = prefilter
    )

    # Prefilter very low quailty cells, if desired
    if (isTRUE(prefilter)) {
        counts <- counts[, rownames(colData), drop = FALSE]
    }

    # Bind the `nCount` column into the colData
    nCount <- cbData[rownames(colData), "nCount", drop = FALSE]
    colData <- cbind(nCount, colData)

    # Cell to sample mappings ==================================================
    cell2sample <- mapCellsToSamples(
        cells = rownames(colData),
        samples = rownames(sampleData)
    )

    # Metadata =================================================================
    metadata <- list(
        "version" = packageVersion,
        "pipeline" = pipeline,
        "level" = level,
        "uploadDir" = uploadDir,
        "sampleDirs" = sampleDirs,
        "sampleMetadataFile" = as.character(sampleMetadataFile),
        "interestingGroups" = interestingGroups,
        "organism" = organism,
        "genomeBuild" = as.character(genomeBuild),
        "ensemblRelease" = as.integer(ensemblRelease),
        "rowRangesMetadata" = rowRangesMetadata,
        "sampleData" = sampleData,
        "cell2sample" = cell2sample,
        "umiType" = umiType,
        "allSamples" = allSamples,
        "prefilter" = prefilter,
        # bcbio pipeline-specific ==============================================
        "projectDir" = projectDir,
        "template" = template,
        "runDate" = runDate,
        "yaml" = yaml,
        "gffFile" = as.character(gffFile),
        "tx2gene" = tx2gene,
        "dataVersions" = dataVersions,
        "programVersions" = programVersions,
        "bcbioLog" = bcbioLog,
        "bcbioCommandsLog" = bcbioCommandsLog,
        "cellularBarcodes" = cbList,
        "cellularBarcodeCutoff" = cellularBarcodeCutoff,
        "loadSingleCell" = match.call()
    )
    # Add user-defined custom metadata, if specified
    if (length(dots)) {
        assert_are_disjoint_sets(names(metadata), names(dots))
        metadata <- c(metadata, dots)
    }

    # Return ===================================================================
    .new.bcbioSingleCell(
        assays = list("raw" = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        isSpike = isSpike
    )
}
