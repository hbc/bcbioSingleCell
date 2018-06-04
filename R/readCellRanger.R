#' Read 10X Genomics Cell Ranger Data
#'
#' Read [10x Genomics Chromium](https://www.10xgenomics.com/software/) cell
#' counts from `barcodes.tsv`, `genes.tsv`, and `matrix.mtx` files.
#'
#' @details This function is a variant of the main [bcbioSingleCell()]
#'   constructor, but optimized for handling Cell Ranger output.
#'
#' @section Directory Structure:
#' Cell Ranger can vary in its output directory structure, but we're requiring a
#' single, consistent data structure for datasets containing multiple samples.
#' Note that Cell Ranger data may not always contain per sample subdirectories,
#' or the "outs" subdirectory. We may make this more flexible in the future, but
#' for now we're making this strict to ensure reproducibility.
#'
#' \preformatted{
#' file.path(
#'     "<uploadDir>",
#'     "<sampleName>",
#'     "outs",
#'     "filtered_gene_bc_matrices*",
#'     "outs",
#'     "<genomeBuild>",
#'     "matrix.mtx"
#' )
#' }
#'
#' @section Sample metadata:
#' A user-defined sample metadata file (`sampleMetadataFile`) is required for
#' multiplexed datasets. Otherwise this can be left `NULL`, and minimal sample
#' data will be used, based on the directory names.
#'
#' @section Reference Data:
#' We strongly recommend supplying the corresponding reference data required for
#' Cell Ranger with the `refdataDir` argument. When set, the function will
#' detect the `organism`, `ensemblRelease`, and `genomeBuild` automatically,
#' based on the 10X `refdataDir` YAML metadata. Additionally, it will convert
#' the gene annotations defined in the GTF file into `GRanges`, which get
#' slotted in [rowRanges()]. Otherwise, the function will attempt to use the
#' most current annotations available from Ensembl, and some gene IDs may not
#' match, due to deprecation in the current Ensembl release.
#'
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams bcbioSingleCell
#' @inheritParams general
#' @param uploadDir Path to Cell Ranger output directory. This directory path
#'   must contain `filtered_gene_bc_matrices*` as a child directory.
#' @param refdataDir Directory path to Cell Ranger reference annotation data.
#'
#' @return `SingleCellExperiment`.
#' @export
#'
#' @examples
#' uploadDir <- system.file("extdata/cellranger", package = "bcbioSingleCell")
#' x <- readCellRanger(uploadDir)
#' show(x)
readCellRanger <- function(
    uploadDir,
    sampleMetadataFile = NULL,
    refdataDir = NULL,
    interestingGroups = "sampleName",
    transgeneNames = NULL,
    spikeNames = NULL,
    ...
) {
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    uploadDir <- normalizePath(uploadDir, winslash = "/", mustWork = TRUE)
    assertIsAStringOrNULL(sampleMetadataFile)
    assertIsAStringOrNULL(refdataDir)
    if (is_a_string(refdataDir)) {
        assert_all_are_dirs(refdataDir)
        refdataDir <- normalizePath(refdataDir, winslash = "/", mustWork = TRUE)
    }
    assert_is_character(interestingGroups)
    assert_is_any_of(transgeneNames, c("character", "NULL"))
    assert_is_any_of(spikeNames, c("character", "NULL"))
    dots <- list(...)
    pipeline <- "cellranger"
    level <- "genes"
    umiType <- "chromium"

    # Legacy arguments =========================================================
    # annotable
    if ("annotable" %in% names(call)) {
        stop("`annotable` argument is defunct")
    }
    dots <- Filter(Negate(is.null), dots)

    # Directory paths ==========================================================
    matrixFiles <- list.files(
        path = uploadDir,
        pattern = "matrix.mtx",
        include.dirs = FALSE,
        full.names = TRUE,
        recursive = TRUE
    )
    # Subset to only include `filtered_gene_bc_matrices*`. Note that
    # aggregation output is labeled `filtered_gene_bc_matrices_mex` by
    # default.
    matrixFiles <- matrixFiles[grepl("filtered_gene_bc_matrices", matrixFiles)]
    assert_is_non_empty(matrixFiles)

    # Combined or per sample matrix support
    if (is_a_string(matrixFiles)) {
        sampleDirs <- uploadDir
    } else {
        message(paste(length(matrixFiles), "samples detected"))
        message(paste(
            "Required per sample output directory structure:",
            file.path(
                "<uploadDir>",
                "<sampleName>",
                "outs",
                "filtered_gene_bc_matrices*",
                "<genomeBuild>",
                "matrix.mtx"
            ),
            sep = "\n"
        ))
        assert_all_are_matching_regex(
            x = matrixFiles,
            pattern = file.path(
                paste0("^", uploadDir),
                "[^/]+", # sampleName
                "outs",
                "filtered_gene_bc_matrices([^/]+)?",
                "[^/]+",
                paste0("matrix.mtx", "$")
            )
        )
        # Sample directories nest the matrix files 4 levels deep
        sampleDirs <- matrixFiles %>%
            dirname() %>%
            dirname() %>%
            dirname() %>%
            dirname()
    }

    names(sampleDirs) <- makeNames(basename(sampleDirs), unique = TRUE)

    # Sample metadata ==========================================================
    if (is_a_string(sampleMetadataFile)) {
        sampleData <- readSampleData(sampleMetadataFile)
    } else {
        sampleData <- minimalSampleData(basename(sampleDirs))
    }

    # Interesting groups =======================================================
    # Ensure internal formatting in camelCase
    interestingGroups <- camel(interestingGroups, strict = FALSE)
    assertFormalInterestingGroups(sampleData, interestingGroups)

    # Subset sample directories by metadata ====================================
    if (nrow(sampleData) < length(sampleDirs)) {
        message("Loading a subset of samples, defined by the metadata file")
        allSamples <- FALSE
        sampleDirs <- sampleDirs %>%
            .[names(sampleDirs) %in% rownames(sampleData)]
        message(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Gene annotations =========================================================
    refJSON <- NULL
    ensemblRelease <- NULL
    rowRangesMetadata <- NULL

    # Stop on multiple genomes (not supported in a single SCE object)
    genomeBuild <- basename(dirname(matrixFiles))
    assert_is_a_string(genomeBuild)
    organism <- detectOrganism(genomeBuild)
    assert_is_a_string(organism)

    # Prepare gene annotations as GRanges
    if (is_a_string(refdataDir)) {
        # JSON data
        refJSONFile <- file.path(refdataDir, "reference.json")
        assert_all_are_existing_files(refJSONFile)
        refJSON <- read_json(refJSONFile)
        # Convert the GTF file to GRanges
        gtfFile <- refJSON[["input_gtf_files"]]
        assert_is_a_string(gtfFile)
        rowRanges <- makeGRangesFromGTF(gtfFile)
        # Get the Ensembl version from the GTF file name.
        # Example: "Homo_sapiens.GRCh37.82.filtered.gtf"
        ensemblRelease <- gtfFile %>%
            str_split("\\.", simplify = TRUE) %>%
            .[1L, 3L] %>%
            as.integer()
    } else {
        # CellRanger uses Ensembl refdata internally. Here we're fetching the
        # annotations with AnnotationHub rather than pulling from the GTF file
        # in the refdata directory. It will also drop genes that are now dead in the
        # current Ensembl release. Don't warn about old Ensembl release version.
        ah <- suppressWarnings(makeGRangesFromEnsembl(
            organism = organism,
            format = level,
            genomeBuild = genomeBuild,
            metadata = TRUE
        ))
        assert_is_list(ah)
        assert_are_identical(names(ah), c("data", "metadata"))
        rowRanges <- ah[["data"]]
        assert_is_all_of(rowRanges, "GRanges")
        rowRangesMetadata <- ah[["metadata"]]
        assert_is_data.frame(rowRangesMetadata)
    }

    # Require gene-to-symbol mappings
    assert_is_subset(
        x = c("geneID", "geneName"),
        y = names(mcols(rowRanges))
    )

    rowData <- as.data.frame(rowRanges)
    rownames(rowData) <- names(rowRanges)

    # Counts ===================================================================
    message("Reading counts at gene level")
    sparseCountsList <- .sparseCountsList(
        sampleDirs = sampleDirs,
        pipeline = pipeline,
        umiType = umiType
    )
    counts <- do.call(cbind, sparseCountsList)

    # Multiplexed sample check =================================================
    # Check to see if multiplexed samples are present and require metadata
    multiplexedPattern <- "^(.+)_(\\d+)_([ACGT]+)$"
    if (any(grepl(multiplexedPattern, colnames(counts)))) {
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

    # Column data ==============================================================
    # Always prefilter, removing very low quality cells with no UMIs or genes
    metrics <- metrics(counts, rowData = rowData, prefilter = TRUE)

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

    # Metadata =================================================================
    metadata <- list(
        "version" = packageVersion,
        "pipeline" = pipeline,
        "level" = level,
        "uploadDir" = uploadDir,
        "sampleDirs" = sampleDirs,
        "sampleMetadataFile" = as.character(sampleMetadataFile),
        "sampleData" = sampleData,
        "interestingGroups" = interestingGroups,
        "cell2sample" = as.factor(cell2sample),
        "organism" = organism,
        "genomeBuild" = as.character(genomeBuild),
        "ensemblRelease" = as.integer(ensemblRelease),
        "rowRangesMetadata" = rowRangesMetadata,
        "umiType" = umiType,
        "allSamples" = allSamples,
        # cellranger pipeline-specific -----------------------------------------
        "refdataDir" = refdataDir,
        "refJSON" = refJSON,
        "call" = match.call()
    )
    # Add user-defined custom metadata, if specified
    if (length(dots)) {
        assert_are_disjoint_sets(names(metadata), names(dots))
        metadata <- c(metadata, dots)
    }

    # Return ===================================================================
    .new.SingleCellExperiment(
        assays = list("counts" = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
}
