#' Coerce Object
#'
#' @rdname coerce
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @param from Class for which the coerce method will perform coercion.
#'
#' @seealso `help(topic = "coerce", package = "methods")`.
#'
#' @examples
#' \dontrun{
#' seurat <- as(bcb, "seurat")
#' }
NULL



# Constructors ====
#' Coerce Legacy bcbio Object to `bcbioSingleCell` class
#'
#' Compatible with old versions created by bcbioSingleCell package.
#' The previous bcbioSinglecell package (note lowercase "c") must be reinstalled
#' to load objects from versions <= 0.0.16.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom utils packageVersion
#'
#' @return [bcbioSingleCell] object.
.coerceLegacy <- function(from) {
    # Check for version
    version <- metadata(from)[["version"]]
    if (is.null(version)) {
        stop(paste(
            "Unknown bcbio object version.",
            "Please reload with 'loadSingleCell()'."
        ), call. = FALSE)
    }
    message(paste(
        paste("Upgrading from", version, "to",
              packageVersion("bcbioSingleCell")),
        paste("Existing metadata:", toString(names(metadata(from)))),
        sep = "\n"
    ))

    se <- as(from, "SummarizedExperiment")
    to <- new("bcbioSingleCell", se)
    bcbio(to) <- bcbio(from)
    validObject(to)

    # Recalculate the cellular barcode metrics
    colData(to) <- calculateMetrics(assay(to), annotable = annotable(to))

    # Update the automatic metadata slots
    metadata(to)[["version"]] <- packageVersion("bcbioSingleCell")
    metadata(to)[["originalVersion"]] <- metadata(from)[["version"]]
    metadata(to)[["upgradeDate"]] <- Sys.Date()

    # Version-specific modifications
    if (version <= package_version("0.0.18")) {
        # Remove GTF file, if present (too large)
        metadata(to)[["gtf"]] <- NULL
    }

    to
}



#' Coerce `bcbioSingleCell` to `seurat`
#'
#' Last tested against CRAN version 2.0.1
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom Seurat CreateSeuratObject FindVariableGenes NormalizeData
#'   ScaleData
#'
#' @return [seurat] object.
.coerceToSeurat <- function(from) {
    # Require that technical replicates are aggregated
    if (any(grepl(x = metadata(from)[["sampleMetadata"]][["sampleID"]],
                  pattern = lanePattern)) |
        "sampleNameAggregate" %in%
        colnames(metadata(from)[["sampleMetadata"]])) {
        stop(paste("'aggregateReplicates()' required",
                   "to merge technical replicates prior to seurat coercion"
        ), call. = FALSE)
    }

    # Require filtered cells and genes only
    from <- .applyFilterCutoffs(from)
    filterParams <- metadata(from)[["filterParams"]]
    print(filterParams)

    # Create the initial `seurat` object
    counts <- counts(from, gene2symbol = TRUE)
    metrics <- metrics(from)

    seurat <- CreateSeuratObject(
        raw.data = counts,
        project = "bcbioSingleCell",
        # Already applied filtering cutoffs for cells and genes
        min.cells = 0,
        min.genes = 0,
        # Default for UMI datasets
        is.expr = 0,
        meta.data = metrics)

    # Check that the dimensions match exactly
    if (!identical(dim(from), dim(slot(seurat, "raw.data")))) {
        stop("Dimension mismatch between bcbioSingleCell and seurat objects",
             call. = FALSE)
    }

    # Stash bcbio run metadata into `misc` slot
    slot(seurat, "misc")[["bcbio"]] <- metadata(from)

    # Normalize and scale the seurat object
    seurat <- seurat %>%
        NormalizeData(
            object = .,
            normalization.method = "LogNormalize",
            scale.factor = 10000) %>%
        FindVariableGenes(
            object = .,
            mean.function = ExpMean,
            dispersion.function = LogVMR,
            do.plot = FALSE) %>%
        ScaleData(
            object = .,
            model.use = "linear")

    print(seurat)
    seurat
}



.coerceToSummarizedExperiment <- function(from) {
    to <- new("SummarizedExperiment")
    slot(to, "colData") <- slot(from, "colData")
    slot(to, "assays") <- slot(from, "assays")
    slot(to, "NAMES") <- slot(from, "NAMES")
    slot(to, "elementMetadata") <- slot(from, "elementMetadata")
    slot(to, "metadata") <- slot(from, "metadata")
    validObject(to)
    to
}



# Methods ====
#' @rdname coerce
#' @name upgrade-bcbioSingleCell
#' @section Upgrade [bcbioSingleCell] to current version:
#' This method adds support for upgrading `bcbioSCDataSet` objects to the latest
#' [bcbioSingleCell] class version. This should be backwards compatible to
#' [bcbioSingleCell] version 0.0.17. Previous objects saved using
#' `bcbioSinglecell` (note case) will likely fail to load with newer versions of
#' the package.
setAs(
    from = "bcbioSCDataSet",
    to = "bcbioSingleCell",
    .coerceLegacy)



#' @rdname coerce
#' @name coerce-bcbioSingleCell-seurat
#' @section [bcbioSingleCell] to [seurat]:
#' Interally, this begins by calling [Seurat::CreateSeuratObject()] without any
#' additional filtering cutoffs, since we already applied them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new [seurat] class object, using [as()] object
#' coercion. Next, global-scaling normalization is applied to the raw counts
#' with [Seurat::NormalizeData()], which (1) normalizes the gene expression
#' measurements for each cell by the total expression, (2) multiplies this by a
#' scale factor (10,000 by default), and (3) log-transforms the result.
#' [Seurat::FindVariableGenes()] is then called, which calculates the average
#' expression and dispersion for each gene, places these genes into bins, and
#' then calculates a z-score for dispersion within each bin. This helps control
#' for the relationship between variability and average expression. Finally, the
#' genes are scaled and centered using the [Seurat::ScaleData()] function.
setAs(
    from = "bcbioSingleCell",
    to = "seurat",
    .coerceToSeurat)



#' @rdname coerce
#' @name coerce-bcbioSingleCell-SummarizedExperiment
#' @section [bcbioSingleCell] to [SummarizedExperiment]:
#' Since [bcbioSingleCell] is an extension of [SummarizedExperiment], this
#' coercion method is very simple. Here we're simply dropping our `@bcbio` slot,
#' which contains raw cellular barcodes and other bcbio-specific metadata.
setAs(
    from = "bcbioSingleCell",
    to = "SummarizedExperiment",
    .coerceToSummarizedExperiment)
