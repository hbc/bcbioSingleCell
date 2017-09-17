#' Load 10X Genomics CellRanger Data
#'
#' Read [10x Genomics Chromium](https://www.10xgenomics.com/software/) cell
#' counts from `barcodes.tsv`, `genes.tsv`, and `matrix.mtx` files.
#'
#' @details This function is a simplified version of [loadSingleCellRun()]
#'   optimized for handling CellRanger output.
#'
#' @author Michael Steinbaugh
#'
#' @inheritParams loadSingleCellRun
#' @param uploadDir Path to CellRanger output directory. This directory path
#'   must contain `filtered_gene_bc_matrices/` as a child.
#' @param refDataDir Directory path to cellranger reference annotation data.
#'
#' @return [bcbioSCDataSet].
#' @export
loadCellRanger <- function(
    uploadDir,
    refDataDir,
    sampleMetadataFile,
    interestingGroups = "sampleName",
    ...) {
    # Initial run setup ====
    pipeline <- "cellranger"
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

    # Reference data ====
    # JSON data
    refDataDir <- normalizePath(refDataDir)
    refJSONFile <- file.path(refDataDir, "reference.json")
    if (!file.exists(refJSONFile)) {
        stop("reference.json file missing")
    }
    refJSON <- read_json(refJSONFile)
    genomeBuild <- refJSON %>%
        .[["genomes"]] %>%
        .[[1L]]
    organism <- detectOrganism(genomeBuild)

    # GTF
    gtfFile <- file.path(refDataDir, "genes", "genes.gtf")
    if (!file.exists(gtfFile)) {
        stop("GTF file missing")
    }
    gtf <- readGTF(gtfFile)

    # gene2symbol mappings
    gene2symbol <- gene2symbolFromGTF(gtf)

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
    metadata <- SimpleList(
        version = packageVersion("bcbioSingleCell"),
        pipeline = pipeline,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = sampleMetadataFile,
        sampleMetadata = sampleMetadata,
        interestingGroups = interestingGroups,
        genomeBuild = genomeBuild,
        organism = organism,
        annotable = annotable,
        gtfFile = gtfFile,
        gtf = gtf,
        gene2symbol = gene2symbol,
        umiType = "chromium",
        allSamples = allSamples,
        multiplexedFASTQ = FALSE,
        # cellranger pipeline-specific
        refDataDir = refDataDir,
        refJSON = refJSON)
    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0L) {
        metadata <- c(metadata, dots)
    }

    # bcbioSCDataSet ===========================================================
    se <- prepareSummarizedExperiment(
        sparseCounts,
        colData = metrics,
        rowData = annotable,
        metadata = metadata)
    new("bcbioSCDataSet", se)
}
