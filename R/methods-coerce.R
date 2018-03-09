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
#' load(system.file("extdata/filtered.rda", package = "bcbioSingleCell"))
#'
#' # Coerce bcbioSingleCell to seurat
#' as(filtered, "seurat")
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
#' @importFrom Seurat CreateSeuratObject
#'
#' @return `seurat`.
.coerceToSeurat <- function(from) {
    # Require that technical replicates are aggregated
    if ("sampleNameAggregate" %in% colnames(sampleData(from))) {
        abort(paste(
            "`aggregateReplicates()` required",
            "to merge technical replicates prior to seurat coercion"
        ))
    }

    # Require filtered cells and genes only
    from <- .applyFilterCutoffs(from)

    # Create the initial `seurat` object
    counts <- counts(from, gene2symbol = TRUE)
    metadata <- colData(from)

    seurat <- CreateSeuratObject(
        raw.data = counts,
        project = "bcbioSingleCell",
        # Already applied filtering cutoffs for cells and genes
        min.cells = 0L,
        min.genes = 0L,
        # Default for UMI datasets
        is.expr = 0L,
        meta.data = metrics
    )

    # Check that the dimensions match exactly
    if (!identical(dim(from), dim(slot(seurat, "raw.data")))) {
        abort("Dimension mismatch between bcbioSingleCell and seurat objects")
    }

    # Stash bcbio run metadata into `misc` slot
    slot(seurat, "misc")[["bcbio"]] <- metadata(from)

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
