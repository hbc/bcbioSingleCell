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
#' @importFrom basejump annotable camel detectOrganism gene2symbolFromGTF
#'   prepareSummarizedExperiment readGTF readSampleMetadataFile
#' @importFrom dplyr mutate
#' @importFrom jsonlite read_json
#' @importFrom magrittr set_colnames
#' @importFrom Matrix cBind
#' @importFrom pbapply pblapply
#' @importFrom stats setNames
#' @importFrom stringr str_split
#'
#' @inherit loadSingleCell
#'
#' @param uploadDir Path to CellRanger output directory. This directory path
#'   must contain `filtered_gene_bc_matrices/` as a child.
#' @param refDataDir Directory path to cellranger reference annotation data.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' extdataDir <- system.file("extdata", package = "bcbioSingleCell")
#' uploadDir <- file.path(extdataDir, "cellranger")
#' refDataDir <- file.path(extdataDir, "refdata-cellranger-hg19-1.2.0")
#' sampleMetadataFile <- file.path(extdataDir, "cellranger.csv")
#' loadCellRanger(
#'     uploadDir = uploadDir,
#'     refDataDir = refDataDir,
#'     sampleMetadataFile = sampleMetadataFile)
#' }
loadCellRanger <- function(
    uploadDir,
    refDataDir,
    interestingGroups = "sampleName",
    sampleMetadataFile = NULL,
    annotable,
    prefilter = TRUE,
    ...) {
    pipeline <- "cellranger"
    umiType <- "chromium"

    # Directory paths ==========================================================
    # Check connection to final upload directory
    if (!dir.exists(uploadDir)) {
        stop("Final upload directory does not exist", call. = FALSE)
    }
    uploadDir <- normalizePath(uploadDir)
    sampleDirs <- .sampleDirs(uploadDir, pipeline = pipeline)

    # Sample metadata ==========================================================
    if (!is.null(sampleMetadataFile)) {
        sampleMetadataFile <- normalizePath(sampleMetadataFile)
        sampleMetadata <- readSampleMetadataFile(sampleMetadataFile)
    } else {
        warning(paste(
            "'sampleMetadataFile' not specified.",
            "Generating minimal sample metadata."
        ), call. = FALSE)
        sampleMetadata <- data.frame(row.names = names(sampleDirs))
        for (i in seq_along(metadataPriorityCols)) {
            sampleMetadata[, metadataPriorityCols[[i]]] <-
                as.factor(names(sampleDirs))
        }
    }
    # Check that `description` matches `sampleDirs`
    if (!all(sampleMetadata[["description"]] %in% basename(sampleDirs))) {
        warning("Sample directory names don't match the sample metadata file",
             call. = FALSE)
    }

    # Interesting groups =======================================================
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    # Check to ensure interesting groups are defined
    if (!all(interestingGroups %in% colnames(sampleMetadata))) {
        stop("Interesting groups missing in sample metadata", call. = FALSE)
    }

    # Subset sample directories by metadata ====================================
    if (nrow(sampleMetadata) < length(sampleDirs)) {
        message("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs %>%
            .[names(sampleDirs) %in% rownames(sampleMetadata)]
        message(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Reference data ===========================================================
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
        .[1, 3] %>%
        as.integer()
    message(paste0("Genome: ", organism, " (", genomeBuild, ")"))
    message(paste("Ensembl release:", ensemblVersion))

    # Gene annotations =========================================================
    if (missing(annotable)) {
        annotable <- basejump::annotable(
            organism,
            genomeBuild = genomeBuild,
            release = ensemblVersion)
    } else if (is.data.frame(annotable)) {
        annotable <- annotable(annotable)
    } else {
        annotable <- NULL
    }
    gtfFile <- file.path(refDataDir, "genes", "genes.gtf")
    if (!file.exists(gtfFile)) {
        # Provide a fallback (for minimal unit testing)
        warning("Reference GTF file missing", call. = FALSE)
        gene2symbol <- annotable(
            organism,
            genomeBuild = genomeBuild,
            format = "gene2symbol")
    } else {
        gtf <- readGTF(gtfFile)
        gene2symbol <- gene2symbolFromGTF(gtf)
    }

    # Counts ===================================================================
    message("Reading counts")
    # Migrate this to `mapply()` method in future update
    sparseList <- pblapply(seq_along(sampleDirs), function(a) {
        .readSparseCounts(
            sampleDirs[a],
            pipeline = pipeline,
            umiType = umiType)
    })
    names(sparseList) <- names(sampleDirs)
    # Cell Ranger outputs at gene-level
    counts <- do.call(Matrix::cBind, sparseList)

    # Metrics ==================================================================
    metrics <- calculateMetrics(
        counts,
        annotable = annotable,
        prefilter = prefilter)
    if (isTRUE(prefilter)) {
        # Subset the counts matrix to match the metrics
        counts <- counts[, rownames(metrics)]
    }

    # Cell to sample mappings ==================================================
    # Check for multiplexed samples. CellRanger outputs these with a trailing
    # number (e.g. `-2$`, which we're sanitizing to `_2$`).
    if (any(grepl(x = colnames(counts), pattern = "_2$"))) {
        if (!"sequence" %in% colnames(sampleMetadata)) {
            stop(paste(
                "'sequence' column must be defined using",
                "'sampleMetadataFile' for multiplexed samples"
            ), call. = FALSE)
        }
        # in the metadata
        map <- str_match(
            string = colnames(counts),
            pattern = "^(.+)_([ACGT]+)_(\\d+)$"
        ) %>%
            as.data.frame() %>%
            set_colnames(
                c("cellID",
                  "description",
                  "barcode",
                  "sequence")) %>%
            mutate_all(as.factor) %>%
            # Note that we can't use minimal sample metadata here
            left_join(sampleMetadata, by = c("description", "sequence"))
        cell2sample <- map[["sampleID"]]
        names(cell2sample) <- map[["cellID"]]
    } else {
    cell2sample <- .cell2sample(
        cells = colnames(counts),
        samples = rownames(sampleMetadata))
    }

    # Metadata =================================================================
    metadata <- list(
        version = packageVersion,
        pipeline = pipeline,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = sampleMetadataFile,
        sampleMetadata = sampleMetadata,
        interestingGroups = interestingGroups,
        cell2sample = cell2sample,
        organism = organism,
        genomeBuild = genomeBuild,
        ensemblVersion = ensemblVersion,
        annotable = annotable,
        gtfFile = gtfFile,
        gene2symbol = gene2symbol,
        umiType = umiType,
        allSamples = allSamples,
        prefilter = prefilter,
        # cellranger pipeline-specific
        refDataDir = refDataDir,
        refJSON = refJSON)
    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0) {
        metadata <- c(metadata, dots)
    }

    # SummarizedExperiment =====================================================
    se <- prepareSummarizedExperiment(
        assays = list(assay = counts),
        rowData = annotable,
        colData = metrics,
        metadata = metadata)
    new("bcbioSingleCell", se)
}
