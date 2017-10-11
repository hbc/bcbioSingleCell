#' Load 10X Genomics CellRanger Data
#'
#' Read [10x Genomics Chromium](https://www.10xgenomics.com/software/) cell
#' counts from `barcodes.tsv`, `genes.tsv`, and `matrix.mtx` files.
#'
#' @details This function is a simplified version of [loadSingleCell()]
#'   optimized for handling CellRanger output.
#'
#' @author Michael Steinbaugh
#'
#' @inherit loadSingleCell
#' @param uploadDir Path to CellRanger output directory. This directory path
#'   must contain `filtered_gene_bc_matrices/` as a child.
#' @param refDataDir Directory path to cellranger reference annotation data.
#'
#' @export
loadCellRanger <- function(
    uploadDir,
    refDataDir,
    interestingGroups = "sampleName",
    sampleMetadataFile,
    prefilter = TRUE,
    ...) {
    pipeline <- "cellranger"
    umiType <- "chromium"
    multiplexedFASTQ <- FALSE

    # Directory paths ====
    # Check connection to final upload directory
    if (!dir.exists(uploadDir)) {
        stop("Final upload directory does not exist", call. = FALSE)
    }
    uploadDir <- normalizePath(uploadDir)
    sampleDirs <- .sampleDirs(uploadDir, pipeline = pipeline)

    # Sample metadata ====
    sampleMetadataFile <- normalizePath(sampleMetadataFile)
    sampleMetadata <- readSampleMetadataFile(sampleMetadataFile)
    # Check that `sampleID` matches `sampleDirs`
    if (!all(sampleMetadata[["sampleID"]] %in% names(sampleDirs))) {
        stop("Sample directory names don't match the sample metadata file",
             call. = FALSE)
    }

    # Interesting groups ====
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    # Check to ensure interesting groups are defined
    if (!all(interestingGroups %in% colnames(sampleMetadata))) {
        stop("Interesting groups missing in sample metadata", call. = FALSE)
    }

    # Subset sample directories by metadata ====
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
        stop("'reference.json' file missing", call. = FALSE)
    }
    refJSON <- read_json(refJSONFile)
    genomeBuild <- refJSON %>%
        .[["genomes"]] %>%
        .[[1]]
    organism <- detectOrganism(genomeBuild)
    # Get the Ensembl version from the JSON reference file (via GTF)
    ensemblVersion <- refJSON %>%
        .[["input_gtf_files"]] %>%
        .[[1]] %>%
        str_split("\\.", simplify = TRUE) %>%
        .[1, 3]
    message(paste("Organism:", organism))
    message(paste("Genome build:", genomeBuild))
    message(paste("Ensembl release:", ensemblVersion))

    # Cell Ranger uses reference GTF file
    gtfFile <- file.path(refDataDir, "genes", "genes.gtf")
    if (!file.exists(gtfFile)) {
        stop("Reference GTF file missing", call. = FALSE)
    }
    gtf <- readGTF(gtfFile)

    # gene2symbol mappings
    gene2symbol <- gene2symbolFromGTF(gtf)

    # Row data =================================================================
    annotable <- annotable(organism, release = ensemblVersion)

    # Assays ===================================================================
    message("Reading counts")
    # Migrate this to `mapply()` method in future update
    sparseList <- pblapply(seq_along(sampleDirs), function(a) {
        .readSparseCounts(sampleDirs[a], pipeline = pipeline)
    }) %>%
        setNames(names(sampleDirs))
    # Cell Ranger outputs at gene-level
    counts <- do.call(Matrix::cBind, sparseList)

    # Column data ==============================================================
    # Calculate the cellular barcode metrics
    metrics <- calculateMetrics(
        counts,
        annotable = annotable,
        prefilter = prefilter)
    if (isTRUE(prefilter)) {
        # Subset the counts matrix to match the metrics
        counts <- counts[, rownames(metrics)]
    }

    # Metadata =================================================================
    metadata <- list(
        version = packageVersion("bcbioSingleCell"),
        pipeline = pipeline,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = sampleMetadataFile,
        sampleMetadata = sampleMetadata,
        interestingGroups = interestingGroups,
        organism = organism,
        genomeBuild = genomeBuild,
        ensemblVersion = ensemblVersion,
        annotable = annotable,
        gtfFile = gtfFile,
        gene2symbol = gene2symbol,
        umiType = umiType,
        allSamples = allSamples,
        multiplexedFASTQ = multiplexedFASTQ,
        prefilter = prefilter,
        # cellranger pipeline-specific
        refDataDir = refDataDir,
        refJSON = refJSON)
    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0) {
        metadata <- c(metadata, dots)
    }

    # Return `bcbioSingleCell` object ==========================================
    se <- prepareSummarizedExperiment(
        assays = list(assay = counts),
        rowData = annotable,
        colData = metrics,
        metadata = metadata)
    new("bcbioSingleCell", se)
}
