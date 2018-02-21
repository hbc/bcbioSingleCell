#' Load 10X Genomics CellRanger Data
#'
#' Read [10x Genomics Chromium](https://www.10xgenomics.com/software/) cell
#' counts from `barcodes.tsv`, `genes.tsv`, and `matrix.mtx` files.
#'
#' @details This function is a simplified version of [loadSingleCell()]
#'   optimized for handling CellRanger output.
#'
#' @note Unlike [loadSingleCell()], the `organism`, `ensemblVersion`, and
#'   `genomeBuild` are always detected automatically, based on the `refdataDir`
#'   YAML metadata. Therefore, these parameters cannot be set by the user.
#'
#' @author Michael Steinbaugh
#'
#' @importFrom basejump annotable camel detectOrganism gene2symbol
#'   gene2symbolFromGTF readGTF
#' @importFrom bcbioBase prepareSummarizedExperiment readSampleMetadataFile
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
#'   must contain `filtered_gene_bc_matrices*` as a child directory.
#' @param refdataDir Directory path to cellranger reference annotation data.
#'
#' @return [bcbioSingleCell].
#' @export
#'
#' @examples
#' extdataDir <- system.file("extdata", package = "bcbioSingleCell")
#' uploadDir <- file.path(extdataDir, "cellranger")
#' refdataDir <- file.path(extdataDir, "refdata-cellranger-hg19-1.2.0")
#' sampleMetadataFile <- file.path(extdataDir, "cellranger.csv")
#' bcb <- loadCellRanger(
#'     uploadDir = uploadDir,
#'     refdataDir = refdataDir,
#'     sampleMetadataFile = sampleMetadataFile)
#' print(bcb)
loadCellRanger <- function(
    uploadDir,
    refdataDir,
    interestingGroups = "sampleName",
    sampleMetadataFile = NULL,
    annotable = TRUE,
    prefilter = TRUE,
    ...) {
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    uploadDir <- normalizePath(uploadDir)
    assert_is_a_string(refdataDir)
    assert_all_are_dirs(refdataDir)
    refdataDir <- normalizePath(refdataDir)
    assert_is_character(interestingGroups)
    assertIsAStringOrNULL(sampleMetadataFile)
    assert_is_any_of(annotable, c("data.frame", "logical", "NULL"))
    if (is.data.frame(annotable)) {
        assertIsAnnotable(annotable)
    }
    assert_is_a_bool(prefilter)

    pipeline <- "cellranger"
    umiType <- "chromium"

    # Directory paths ==========================================================
    sampleDirs <- .sampleDirs(uploadDir, pipeline = pipeline)

    # Sample metadata ==========================================================
    if (is_a_string(sampleMetadataFile)) {
        sampleMetadata <- readSampleMetadataFile(sampleMetadataFile)
    } else {
        warn(paste(
            "`sampleMetadataFile` not specified.",
            "Generating minimal sample metadata."
        ))
        sampleMetadata <- data.frame(row.names = names(sampleDirs))
        for (i in seq_along(metadataPriorityCols)) {
            sampleMetadata[, metadataPriorityCols[[i]]] <-
                as.factor(names(sampleDirs))
        }
    }
    assert_is_subset(sampleMetadata[["description"]], basename(sampleDirs))

    # Interesting groups =======================================================
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    assertFormalInterestingGroups(sampleMetadata, interestingGroups)

    # Subset sample directories by metadata ====================================
    if (nrow(sampleMetadata) < length(sampleDirs)) {
        inform("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs %>%
            .[names(sampleDirs) %in% rownames(sampleMetadata)]
        inform(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Reference data ===========================================================
    # JSON data
    refJSONFile <- file.path(refdataDir, "reference.json")
    assert_all_are_existing_files(refJSONFile)
    refJSON <- read_json(refJSONFile)
    genomeBuild <- refJSON %>%
        .[["genomes"]] %>%
        .[[1L]]
    # Only a single genome is currently supported
    assert_is_a_string(genomeBuild)
    organism <- detectOrganism(genomeBuild)
    assert_is_a_string(organism)
    # Get the Ensembl version from the JSON reference file (via GTF)
    ensemblVersion <- refJSON %>%
        .[["input_gtf_files"]] %>%
        .[[1L]] %>%
        str_split("\\.", simplify = TRUE) %>%
        .[1L, 3L] %>%
        as.integer()
    assert_is_an_integer(ensemblVersion)
    assert_all_are_positive(ensemblVersion)
    if (ensemblVersion < 87L) {
        release <- NULL
    } else {
        release <- ensemblVersion
    }
    inform(paste(
        paste("Organism:", organism),
        paste("Genome build:", genomeBuild),
        paste("Ensembl release:", ensemblVersion),
        sep = "\n"
    ))

    # Gene annotations =========================================================
    if (isTRUE(annotable)) {
        annotable <- annotable(
            organism,
            genomeBuild = genomeBuild,
            release = release,
            uniqueSymbol = FALSE)
    } else if (is.data.frame(annotable)) {
        annotable <- annotable(annotable)
    } else {
        warn("Loading run without gene annotations")
        annotable <- NULL
    }

    # GTF annotations
    gtfFile <- file.path(refdataDir, "genes", "genes.gtf")
    if (identical(
        x = refdataDir,
        y = system.file(
            "extdata/refdata-cellranger-hg19-1.2.0",
            package = "bcbioSingleCell")
    )) {
        inform("Minimal working example doesn't contain a GTF file")
        gtf <- NULL
        inform("Obtaining gene2symbol from Ensembl")
        gene2symbol <- gene2symbol(
            organism,
            genomeBuild = genomeBuild,
            release = release)
    } else {
        assert_all_are_existing_files(gtfFile)
        gtf <- readGTF(gtfFile)
        gene2symbol <- gene2symbolFromGTF(gtf)
    }

    # Counts ===================================================================
    inform("Reading counts at gene level")
    sparseCountsList <- .sparseCountsList(
        sampleDirs = sampleDirs,
        pipeline = pipeline,
        umiType = umiType)

    # Cell Ranger always outputs at gene level
    counts <- do.call(Matrix::cBind, sparseCountsList)

    # Metrics ==================================================================
    metrics <- calculateMetrics(
        object = counts,
        annotable = annotable,
        prefilter = prefilter)

    if (isTRUE(prefilter)) {
        # Subset the counts matrix to match the metrics
        counts <- counts[, rownames(metrics), drop = FALSE]
    }

    # Cell to sample mappings ==================================================
    # Check for multiplexed samples. CellRanger outputs these with a trailing
    # number (e.g. `-2$`, which we're sanitizing to `_2$`).
    if (any(grepl(x = colnames(counts), pattern = "_2$"))) {
        if (!"index" %in% colnames(sampleMetadata)) {
            abort(paste(
                "`index` column must be defined using",
                "`sampleMetadataFile` for multiplexed samples"
            ))
        }
        # Prepare data.frame of barcode mappings
        map <- str_match(
            string = colnames(counts),
            pattern = "^(.+)_([ACGT]+)_(\\d+)$"
        ) %>%
            as.data.frame() %>%
            set_colnames(c(
                "cellID",
                "description",
                "barcode",
                "index"
            )) %>%
            mutate_all(as.factor) %>%
            # Note that we can't use minimal sample metadata here
            left_join(sampleMetadata, by = c("description", "index"))
        cell2sample <- map[["sampleID"]]
        names(cell2sample) <- map[["cellID"]]
    } else {
        cell2sample <- mapCellsToSamples(
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
        refdataDir = refdataDir,
        refJSON = refJSON)
    # Add user-defined custom metadata, if specified
    dots <- list(...)
    if (length(dots) > 0L) {
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
