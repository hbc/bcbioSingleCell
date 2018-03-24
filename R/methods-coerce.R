# FIXME Don't use `counts(from, convertGenesToSymbols = TRUE)`


#' Coerce Object
#'
#' @name coerce
#' @author Michael Steinbaugh
#'
#' @param from Class for which the coerce method will perform coercion.
#'
#' @seealso `help(topic = "coerce", package = "methods")`.
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell to seurat
#' as(bcb_small, "seurat")
NULL



# Constructors =================================================================
#' Coerce bcbioSingleCell to seurat
#'
#' Last tested against CRAN version 2.0.1
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @return `seurat`.
.coerceToSeurat <- function(from) {
    # Require that technical replicates are aggregated
    if ("sampleNameAggregate" %in% colnames(sampleData(from))) {
        inform("Use `aggregateReplicates()` to combine technical replicates")
    }

    # Require filtered cells and genes only
    from <- .applyFilterCutoffs(from)

    # Create seurat object
    raw <- counts(from)
    rownames(raw) <- makeNames(gene2symbol(from)[["geneName"]])
    colData <- colData(from, return = "data.frame")
    seurat <- CreateSeuratObject(
        raw.data = raw,
        project = "bcbioSingleCell",
        # Already applied filtering cutoffs for cells and genes
        min.cells = 0L,
        min.genes = 0L,
        # Default for UMI datasets
        is.expr = 0L,
        meta.data = colData
    )

    # Check that the dimensions match exactly
    if (!identical(dim(from), dim(slot(seurat, "raw.data")))) {
        abort("Dimension mismatch between bcbioSingleCell and seurat")
    }

    # Stash bcbio run metadata into `misc` slot
    bcbio <- list(
        "rowRanges" = rowRanges(from),
        "metadata" = metadata(from)
    )
    slot(seurat, "misc")[["bcbio"]] <- bcbio

    seurat
}



# Methods ======================================================================
#' @rdname coerce
#' @name coerce-bcbioSingleCell-seurat
#' @section bcbioSingleCell to seurat:
#' Interally [Seurat::CreateSeuratObject()] is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new `seurat` class object, using [as()] object
#' coercion.
setAs(
    from = "bcbioSingleCell",
    to = "seurat",
    .coerceToSeurat
)
