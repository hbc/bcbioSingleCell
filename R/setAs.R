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
#' load(system.file(
#'     file.path("extdata", "filtered.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # Coerce bcbioSingleCell to seurat
#' as(filtered, "seurat")
NULL



# Constructors =================================================================
#' Coerce Legacy bcbio Object to bcbioSingleCell class
#'
#' Compatible with old versions created by bcbioSingleCell package.
#' The previous bcbioSinglecell package (note lowercase "c") must be reinstalled
#' to load objects from versions <= 0.0.16.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @return [bcbioSingleCell] object.
.coerceLegacy <- function(from) {
    # Check for version
    version <- metadata(from)[["version"]]
    if (is.null(version)) {
        abort(paste(
            "Unknown bcbio object version.",
            "Please reload with 'loadSingleCell()'."
        ))
    }
    inform(paste(
        paste("Upgrading from", version, "to", packageVersion),
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
    metadata(to)[["version"]] <- packageVersion
    metadata(to)[["originalVersion"]] <- metadata(from)[["version"]]
    metadata(to)[["upgradeDate"]] <- Sys.Date()

    # Version-specific modifications
    if (version <= package_version("0.0.18")) {
        # Remove GTF file, if present (too large)
        metadata(to)[["gtf"]] <- NULL
    }

    to
}



#' Coerce bcbioSingleCell to seurat
#'
#' Last tested against CRAN version 2.0.1
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom Seurat CreateSeuratObject
#'
#' @return [seurat] object.
.coerceToSeurat <- function(from) {
    # Require that technical replicates are aggregated
    if ("sampleNameAggregate" %in% colnames(sampleMetadata(from))) {
        abort(paste("`aggregateReplicates()` required",
                   "to merge technical replicates prior to seurat coercion"
        ))
    }

    # Require filtered cells and genes only
    from <- .applyFilterCutoffs(from)

    # Create the initial `seurat` object
    counts <- counts(from, gene2symbol = TRUE)
    metrics <- metrics(from)

    seurat <- CreateSeuratObject(
        raw.data = counts,
        project = "bcbioSingleCell",
        # Already applied filtering cutoffs for cells and genes
        min.cells = 0L,
        min.genes = 0L,
        # Default for UMI datasets
        is.expr = 0L,
        meta.data = metrics)

    # Check that the dimensions match exactly
    if (!identical(dim(from), dim(slot(seurat, "raw.data")))) {
        abort("Dimension mismatch between bcbioSingleCell and seurat objects")
    }

    # Stash bcbio run metadata into `misc` slot
    slot(seurat, "misc")[["bcbio"]] <- metadata(from)

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



# Methods ======================================================================
#' @rdname coerce
#' @name upgrade-bcbioSingleCell
#' @section Upgrade bcbioSingleCell to current version:
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
#' @section bcbioSingleCell to seurat:
#' Interally [Seurat::CreateSeuratObject()] is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new [seurat] class object, using [as()] object
#' coercion.
setAs(
    from = "bcbioSingleCell",
    to = "seurat",
    .coerceToSeurat)



#' @rdname coerce
#' @name coerce-bcbioSingleCell-SummarizedExperiment
#' @section bcbioSingleCell to SummarizedExperiment:
#' Since [bcbioSingleCell] is an extension of [SummarizedExperiment], this
#' coercion method is very simple. Here we're simply dropping our `@bcbio` slot,
#' which contains raw cellular barcodes and other bcbio-specific metadata.
setAs(
    from = "bcbioSingleCell",
    to = "SummarizedExperiment",
    .coerceToSummarizedExperiment)
