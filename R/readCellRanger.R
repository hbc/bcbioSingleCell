#' Read 10X Genomics Cell Ranger Data
#'
#' Read [10x Genomics Chromium](https://www.10xgenomics.com/software/) cell
#' counts from `barcodes.tsv`, `genes.tsv`, and `matrix.mtx` files.
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
#'     "outs",
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
#' @family Read Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams bcbioSingleCell
#' @inheritParams general
#' @param uploadDir Path to Cell Ranger output directory. This directory path
#'   must contain `filtered_gene_bc_matrices*` as a child directory.
#' @param filtered Use filtered (recommended) or raw counts. Note that raw
#'   counts still contain only whitelisted cellular barcodes.
#' @param format Output format, either MatrixMarket ("`mtx`") or HDF5
#'   ("`hdf5`").
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
    format = c("mtx", "hdf5"),
    filtered = TRUE,
    organism = NULL,
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
    format <- match.arg(format)
    assert_is_a_bool(filtered)
    assertIsAStringOrNULL(organism)
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
    stopifnot(!"annotable" %in% names(call))
    dots <- Filter(Negate(is.null), dots)

    # Sample directories =======================================================
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
        sampleDirs <- sampleDirs[rownames(sampleData)]
        message(paste(length(sampleDirs), "samples matched by metadata"))
    } else {
        allSamples <- TRUE
    }

    # Counts ===================================================================
    # This step can be slow over sshfs, recommend running on an HPC
    message("Reading counts at gene level")
    counts <- .readCounts(
        sampleDirs = sampleDirs,
        pipeline = pipeline,
        format = format,
        filtered = filtered
    )

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

    # Row data =================================================================
    refJSON <- NULL
    genomeBuild <- NULL
    ensemblRelease <- NULL
    rowRangesMetadata <- NULL

    # Prepare gene annotations as GRanges
    if (is_a_string(refdataDir)) {
        # JSON data
        refJSONFile <- file.path(refdataDir, "reference.json")
        assert_all_are_existing_files(refJSONFile)
        refJSON <- read_json(refJSONFile)
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
    } else if (is_a_string(organism)) {
        # CellRanger uses Ensembl refdata internally. Here we're fetching the
        # annotations with AnnotationHub rather than pulling from the GTF file
        # in the refdata directory. It will also drop genes that are now dead in
        # the current Ensembl release. Don't warn about old Ensembl release
        # version.
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
    } else {
        rowRanges <- emptyRanges(rownames(counts))
    }
    rowData <- as.data.frame(rowRanges)
    rownames(rowData) <- names(rowRanges)

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
        version = packageVersion,
        pipeline = pipeline,
        level = level,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = as.character(sampleMetadataFile),
        sampleData = sampleData,
        interestingGroups = interestingGroups,
        cell2sample = as.factor(cell2sample),
        organism = organism,
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease),
        rowRangesMetadata = rowRangesMetadata,
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

    # Return ===================================================================
    .new.SingleCellExperiment(
        assays = list(counts = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
}
