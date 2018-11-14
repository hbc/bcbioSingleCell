# FIXME Break these back out into separate files.
# FIXME Either make `sampleName` required or strip it from minimal examples.
# TODO Check to see if we can import tx2gene.csv
# FIXME Can we parse the CellRanger `runDate` from the refData YAML?
# FIXME Allow this function to work if the user points at dir containing matrix.
# FIXME Add documentation about simple mode.



#' @inherit CellRanger-class
#' @export
#'
#' @details
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
#' @author Michael Steinbaugh
#' @export
#'
#' @inheritParams basejump::params
#' @inheritParams bcbioSingleCell
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
CellRanger <- function(  # nolint
    uploadDir,
    format = c("mtx", "hdf5"),
    filtered = TRUE,
    sampleMetadataFile = NULL,
    organism = NULL,
    ensemblRelease = NULL,
    genomeBuild = NULL,
    refdataDir = NULL,
    gffFile = NULL,
    transgeneNames = NULL,
    spikeNames = NULL,
    interestingGroups = "sampleName"
) {
    assert_is_a_string(uploadDir)
    assert_all_are_dirs(uploadDir)
    uploadDir <- realpath(uploadDir)
    format <- match.arg(format)
    assert_is_a_bool(filtered)
    assertIsStringOrNULL(sampleMetadataFile)
    assertIsStringOrNULL(organism)
    assertIsAnImplicitIntegerOrNULL(ensemblRelease)
    assertIsStringOrNULL(genomeBuild)
    assertIsStringOrNULL(refdataDir)
    if (is_a_string(refdataDir)) {
        assert_all_are_dirs(refdataDir)
        refdataDir <- realpath(refdataDir)
    }
    assertIsStringOrNULL(gffFile)
    assert_is_any_of(transgeneNames, c("character", "NULL"))
    assert_is_any_of(spikeNames, c("character", "NULL"))
    assert_is_character(interestingGroups)

    pipeline <- "cellranger"
    level <- "genes"
    umiType <- "chromium"

    # Sample files -------------------------------------------------------------
    sampleFiles <- .sampleFiles.cellranger(uploadDir)

    # Sequencing lanes ---------------------------------------------------------
    lanes <- detectLanes(sampleFiles)

    # Sample metadata ----------------------------------------------------------
    allSamples <- TRUE
    sampleData <- NULL

    if (is_a_string(sampleMetadataFile)) {
        sampleData <- readSampleData(sampleMetadataFile)
        # Allow sample selection by with this file.
        if (nrow(sampleData) < length(sampleFiles)) {
            message("Loading a subset of samples, defined by the metadata.")
            allSamples <- FALSE
            sampleFiles <- sampleFiles[rownames(sampleData)]
            message(paste(length(sampleFiles), "samples matched by metadata."))
        }
    }

    # Counts -------------------------------------------------------------------
    counts <- .import.cellranger(sampleFiles)

    # Row data -----------------------------------------------------------------
    refJSON <- NULL

    # Prepare gene annotations as GRanges.
    if (is_a_string(refdataDir)) {
        message("Using 10X Genomics reference data for annotations.")
        message(paste("refdataDir:", refdataDir))
        # JSON data.
        refJSONFile <- file.path(refdataDir, "reference.json")
        assert_all_are_existing_files(refJSONFile)
        refJSON <- import(refJSONFile)
        # Get the genome build from JSON metadata.
        genomeBuild <- unlist(refJSON[["genomes"]])
        assert_is_a_string(genomeBuild)
        # Convert the GTF file to GRanges.
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
        message("Using `makeGRangesFromEnsembl()` for annotations.")
        rowRanges <- makeGRangesFromEnsembl(
            organism = organism,
            level = level,
            genomeBuild = genomeBuild,
            release = ensemblRelease
        )
        if (is.null(genomeBuild)) {
            genomeBuild <- metadata(rowRanges)[["genomeBuild"]]
        }
        if (is.null(ensemblRelease)) {
            ensemblRelease <- metadata(rowRanges)[["ensemblRelease"]]
        }
    } else {
        message("Unknown organism. Skipping annotations.")
        rowRanges <- emptyRanges(rownames(counts))
    }
    assert_is_all_of(rowRanges, "GRanges")

    # Column data --------------------------------------------------------------
    # Automatic sample metadata.
    if (is.null(sampleData)) {
        # Define the grep pattern to use for sample ID extraction.
        pattern <- "^(.+)_[ACGT]+$"
        if (all(grepl(pattern, colnames(counts)))) {
            match <- str_match(
                string = colnames(counts),
                pattern = pattern
            )
            samples <- unique(match[, 2L, drop = TRUE])
        } else if (has_length(sampleFiles, n = 1L)) {
            samples <- names(sampleFiles)
        }
        sampleData <- minimalSampleData(samples)
    }

    # Always prefilter, removing very low quality cells with no UMIs or genes.
    colData <- metrics.matrix(
        object = counts,
        rowRanges = rowRanges,
        prefilter = TRUE
    )

    # Subset the counts to match the prefiltered metrics.
    assert_is_subset(rownames(colData), colnames(counts))
    counts <- counts[, rownames(colData), drop = FALSE]

    # Join `sampleData` into cell-level `colData`.
    if (has_length(nrow(sampleData), n = 1L)) {
        colData[["sampleID"]] <- as.factor(rownames(sampleData))
    } else {
        colData[["sampleID"]] <- mapCellsToSamples(
            cells = rownames(colData),
            samples = rownames(sampleData)
        )
    }

    # Metadata -----------------------------------------------------------------
    # Interesting groups.
    interestingGroups <- camel(interestingGroups)
    assert_is_subset(interestingGroups, colnames(sampleData))

    metadata <- list(
        version = packageVersion,
        pipeline = pipeline,
        level = level,
        uploadDir = uploadDir,
        sampleDirs = sampleDirs,
        sampleMetadataFile = as.character(sampleMetadataFile),
        interestingGroups = interestingGroups,
        organism = organism,
        genomeBuild = as.character(genomeBuild),
        ensemblRelease = as.integer(ensemblRelease),
        umiType = umiType,
        allSamples = allSamples,
        lanes = lanes,
        # cellranger-specific --------------------------------------------------
        refdataDir = refdataDir,
        refJSON = refJSON,
        call = match.call()
    )

    # Return -------------------------------------------------------------------
    .new.CellRanger(
        assays = list(counts = counts),
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
}



.new.CellRanger <-  # nolint
    function(...) {
        new(
            Class = "CellRanger",
            makeSingleCellExperiment(...)
        )
    }
