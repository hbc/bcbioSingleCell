#' Load 10X Genomics CellRanger Data
#'
#' Read [10x Genomics Chromium](https://www.10xgenomics.com/software/) cell
#' counts from `barcodes.tsv`, `genes.tsv`, and `matrix.mtx` files.
#'
#' @details This function is a simplified version of [loadSingleCellRun()]
#'   optimized for handling CellRanger output.
#'
#' @rdname loadCellRanger
#' @name loadCellRanger
#'
#' @inheritParams loadSingleCellRun
#' @param object Path to CellRanger output directory. This directory path must
#'   contain `filtered_gene_bc_matrices/` as a child.
#'
#' @return [bcbioSCDataSet].
NULL



# Methods ====
#' @rdname loadCellRanger
#' @export
setMethod("loadCellRanger", "character", function(
    object,
    sampleMetadataFile,
    interestingGroups = "sampleName",
    gtfFile = NULL,
    ...) {
    # Initial run setup ====
    pipeline <- "cellranger"
    uploadDir <- object
    if (!dir.exists(uploadDir)) {
        stop("Final upload directory does not exist", call. = FALSE)
    }
    uploadDir <- normalizePath(uploadDir)
    sampleDirs <- .sampleDirs(uploadDir, pipeline = pipeline)

    # Sample metadata ====
    sampleMetadataFile <- normalizePath(sampleMetadataFile)
    sampleMetadata <- .readSampleMetadataFile(
        sampleMetadataFile, sampleDirs, pipeline)

    # Check to ensure interesting groups are defined
    if (!all(interestingGroups %in% colnames(sampleMetadata))) {
        stop("Interesting groups missing in sample metadata")
    }

    # Check to see if a subset of samples is requested via the metadata file.
    # This matches by the reverse complement sequence of the index barcode.
    if (nrow(sampleMetadata) < length(sampleDirs)) {
        message("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs %>%
            .[names(sampleDirs) %in% rownames(sampleMetadata)]
        message(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Get genome build from `sampleDirs`
    genomeBuild <- basename(sampleDirs) %>% unique
    if (length(genomeBuild) > 1L) {
        stop("Multiple genomes detected -- not supported")
    }

    # Row data =================================================================
    annotable <- annotable(genomeBuild)

    # Assays ===================================================================
    message("Reading counts")
    # Migrate this to `mapply()` method in future update
    sparseList <- pblapply(seq_along(sampleDirs), function(a) {
        sparseCounts <- .readSparseCounts(sampleDirs[a], pipeline = pipeline)
        # Pre-filter using cellular barcode summary metrics
        metrics <- calculateMetrics(sparseCounts, annotable)
        sparseCounts[, rownames(metrics)]
    }) %>%
        setNames(names(sampleDirs))
    sparseCounts <- do.call(Matrix::cBind, sparseList)

    # Column data ==============================================================
    metrics <- calculateMetrics(sparseCounts, annotable)

    # Metadata =================================================================
    # gene2symbol mappings ====
    if (!is.null(gtfFile)) {
        gtf <- readGTF(gtfFile)
        gene2symbol <- gene2symbolFromGTF(gtf)
    } else {
        gtf <- NULL
        gene2symbol <- gene2symbol(genomeBuild)
    }

    metadata <- SimpleList(
        version = packageVersion("bcbioSinglecell"),
        pipeline = pipeline,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = sampleMetadataFile,
        sampleMetadata = sampleMetadata,
        interestingGroups = interestingGroups,
        genomeBuild = genomeBuild,
        annotable = annotable,
        gtfFile = gtfFile,
        gtf = gtf,
        gene2symbol = gene2symbol,
        umiType = "chromium",
        allSamples = allSamples,
        multiplexedFASTQ = FALSE)
    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0L) {
        metadata <- c(metadata, dots)
    }

    # bcbioSCDataSet ===========================================================
    se <- prepareSE(
        sparseCounts,
        colData = metrics,
        rowData = annotable,
        metadata = metadata)
    new("bcbioSCDataSet", se)
})
