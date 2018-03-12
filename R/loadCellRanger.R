# TODO Consolidate code where applicable with loadSingleCell
# TODO Build an EnsDb for the old release (84) currently used by cellranger



#' Load 10X Genomics Cell Ranger Data
#'
#' Read [10x Genomics Chromium](https://www.10xgenomics.com/software/) cell
#' counts from `barcodes.tsv`, `genes.tsv`, and `matrix.mtx` files.
#'
#' @details This function is a simplified version of [loadSingleCell()]
#'   optimized for handling Cell Ranger output.
#'
#' @note Unlike [loadSingleCell()], the `organism`, `ensemblRelease`, and
#'   `genomeBuild` are always detected automatically, based on the 10X
#'   `refDataDir` YAML metadata. Therefore, these parameters cannot be set by
#'   the user.
#'
#' @author Michael Steinbaugh
#'
#' @importFrom basejump camel detectOrganism ensembl
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
#' @param uploadDir Path to Cell Ranger output directory. This directory path
#'   must contain `filtered_gene_bc_matrices*` as a child directory.
#' @param refDataDir Directory path to Cell Ranger reference annotation data.
#'
#' @return `bcbioSingleCell`.
#' @export
#'
#' @examples
#' extdataDir <- system.file("extdata", package = "bcbioSingleCell")
#' uploadDir <- file.path(extdataDir, "cellranger")
#' refDataDir <- file.path(extdataDir, "refdata-cellranger-hg19-1.2.0")
#' sampleMetadataFile <- file.path(extdataDir, "cellranger.csv")
#' loadCellRanger(
#'     uploadDir = uploadDir,
#'     refDataDir = refDataDir,
#'     sampleMetadataFile = sampleMetadataFile
#' )
loadCellRanger <- function(
    uploadDir,
    refDataDir,
    interestingGroups = "sampleName",
    sampleMetadataFile = NULL,
    isSpike = NULL,
    prefilter = TRUE,
    ...
) {
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    uploadDir <- normalizePath(uploadDir)
    assert_is_a_string(refDataDir)
    assert_all_are_dirs(refDataDir)
    refDataDir <- normalizePath(refDataDir)
    assert_is_character(interestingGroups)
    assertIsAStringOrNULL(sampleMetadataFile)
    assertIsCharacterOrNULL(isSpike)
    assert_is_a_bool(prefilter)
    dots <- list(...)

    pipeline <- "cellranger"
    level <- "genes"
    umiType <- "chromium"

    # Legacy arguments =========================================================
    # annotable
    if ("annotable" %in% names(call)) {
        abort("User-defined `annotable` no longer allowed")
    }
    # refdataDir
    if ("refdataDir" %in% names(call)) {
        warn("Use `refDataDir` instead of `refdataDir`")
        refDataDir <- call[["refdataDir"]]
        dots[["refdataDir"]] <- NULL
    }
    dots <- Filter(Negate(is.null), dots)

    # Directory paths ==========================================================
    sampleDirs <- .sampleDirs(uploadDir, pipeline = pipeline)

    # Sample metadata ==========================================================
    if (is_a_string(sampleMetadataFile)) {
        sampleData <- readSampleMetadataFile(sampleMetadataFile)
    } else {
        warn(paste(
            "`sampleMetadataFile` not specified.",
            "Generating minimal sample metadata."
        ))
        sampleData <- data.frame(row.names = names(sampleDirs))
        for (i in seq_along(metadataPriorityCols)) {
            sampleData[, metadataPriorityCols[[i]]] <-
                as.factor(names(sampleDirs))
        }
    }
    assert_is_subset(sampleData[["description"]], basename(sampleDirs))

    # Interesting groups =======================================================
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    assertFormalInterestingGroups(sampleData, interestingGroups)

    # Subset sample directories by metadata ====================================
    if (nrow(sampleData) < length(sampleDirs)) {
        inform("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs %>%
            .[names(sampleDirs) %in% rownames(sampleData)]
        inform(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Reference data ===========================================================
    # JSON data
    refJSONFile <- file.path(refDataDir, "reference.json")
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
    ensemblRelease <- refJSON %>%
        .[["input_gtf_files"]] %>%
        .[[1L]] %>%
        str_split("\\.", simplify = TRUE) %>%
        .[1L, 3L] %>%
        as.integer()
    assert_is_an_integer(ensemblRelease)
    assert_all_are_positive(ensemblRelease)
    inform(paste(
        paste("Organism:", organism),
        paste("Genome build:", genomeBuild),
        paste("Ensembl release:", ensemblRelease),
        sep = "\n"
    ))

    # Gene annotations =========================================================
    ahMeta <- NULL
    # ah = AnnotationHub
    ah <- ensembl(
        organism = organism,
        format = level,
        genomeBuild = genomeBuild,
        release = ensemblRelease,
        return = "GRanges",
        metadata = TRUE
    )
    assert_is_list(ah)
    assert_are_identical(names(ah), c("data", "metadata"))
    rowRanges <- ah[["data"]]
    assert_is_all_of(rowRanges, "GRanges")
    ahMeta <- ah[["metadata"]]
    assert_all_are_matching_regex(
        x = ahMeta[["id"]],
        pattern = "^AH\\d+$"
    )

    # GTF annotations
    gtfFile <- file.path(refDataDir, "genes", "genes.gtf")
    if (identical(
        x = refDataDir,
        y = system.file(
            "extdata/refdata-cellranger-hg19-1.2.0",
            package = "bcbioSingleCell"
        )
    )) {
        inform("Minimal working example doesn't contain a GTF file")
        gtf <- NULL
        inform("Obtaining gene2symbol from Ensembl")
        gene2symbol <- gene2symbol(
            organism,
            genomeBuild = genomeBuild,
            release = release
        )
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
        umiType = umiType
    )

    # Cell Ranger always outputs at gene level
    counts <- do.call(Matrix::cBind, sparseCountsList)

    # Metrics ==================================================================
    metrics <- calculateMetrics(
        object = counts,
        rowRanges = rowRanges,
        prefilter = prefilter
    )

    if (isTRUE(prefilter)) {
        # Subset the counts matrix to match the metrics
        counts <- counts[, rownames(metrics), drop = FALSE]
    }

    # Cell to sample mappings ==================================================
    # Check for multiplexed samples. CellRanger outputs these with a trailing
    # number (e.g. `-2$`, which we're sanitizing to `_2$`).
    if (any(grepl(x = colnames(counts), pattern = "_2$"))) {
        if (!"index" %in% colnames(sampleData)) {
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
            left_join(sampleData, by = c("description", "index"))
        cell2sample <- map[["sampleID"]]
        names(cell2sample) <- map[["cellID"]]
    } else {
        cell2sample <- mapCellsToSamples(
            cells = colnames(counts),
            samples = rownames(sampleData)
        )
    }

    # Metadata =================================================================
    metadata <- list(
        "version" = packageVersion,
        "pipeline" = pipeline,
        "uploadDir" = uploadDir,
        "sampleDirs" = sampleDirs,
        "sampleMetadataFile" = as.character(sampleMetadataFile),
        "sampleData" = sampleData,
        "interestingGroups" = interestingGroups,
        "cell2sample" = cell2sample,
        "organism" = organism,
        "genomeBuild" = genomeBuild,
        "ensemblRelease" = ensemblRelease,
        "annotationHub" = as.list(ahMeta),
        "gtfFile" = gtfFile,
        "gene2symbol" = gene2symbol,
        "umiType" = umiType,
        "allSamples" = allSamples,
        "prefilter" = prefilter,
        # cellranger pipeline-specific ====
        "refDataDir" = refDataDir,
        "refJSON" = refJSON,
        "loadCellRanger" = match.call()
    )
    # Add user-defined custom metadata, if specified
    if (length(dots) > 0L) {
        assert_are_disjoint_sets(names(metadata), names(dots))
        metadata <- c(metadata, dots)
    }

    # Return ===================================================================
    # Generate RangedSummarizedExperiment
    rse <- prepareSummarizedExperiment(
        assays = list(assay = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        isSpike = isSpike
    )
    # FIXME `as(rse, "SingleCellExperiment)` coercion doesn't return
    # objectVersion correctly in SingleCellExperiment v1.0.0
    sce <- SingleCellExperiment(
        assays = assays(rse),
        rowRanges = rowRanges(rse),
        colData = colData(rse),
        metadata = metadata(rse)
    )
    # Define spikeNames for spike-in sequences
    if (is.character(isSpike)) {
        for (i in seq_along(isSpike)) {
            isSpike(sce, isSpike[[i]]) <- isSpike[[i]]
        }
    }
    new("bcbioSingleCell", sce)
}
