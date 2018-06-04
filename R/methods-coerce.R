#' Methods for Coercing an Object to a Class
#'
#' @name coerce
#' @aliases as
#' @family S4 Object
#' @author Michael Steinbaugh
#'
#' @inherit bcbioBase::coerce
#'
#' @seealso
#' - [Seurat::CreateSeuratObject()].
#' - [Seurat::Convert()].
#' - [methods::as()].
#' - [methods::coerce()].
NULL



# Methods ======================================================================
#' @rdname coerce
#' @name coerce-SingleCellExperiment-seurat
#'
#' @section SingleCellExperiment to seurat:
#' Interally [Seurat::CreateSeuratObject()] is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new `seurat` class object.
#'
#' @examples
#' # SingleCellExperiment to seurat ====
#' x <- as(indrops_small, "seurat")
#' class(x)
#' print(x)
setAs(
    from = "SingleCellExperiment",
    to = "seurat",
    function(from) {
        # Require that technical replicates are aggregated
        if ("aggregate" %in% colnames(sampleData(from))) {
            message(paste(
                "`aggregate` metadata column detected.",
                "Use `aggregateReplicates()` to combine technical replicates.",
                sep = "\n"
            ))
        }

        # Convert gene identifiers to symbols
        from <- convertGenesToSymbols(from)

        # Create the seurat object
        to <- CreateSeuratObject(
            raw.data = counts(from),
            project = "bcbioSingleCell",
            # Already applied filtering cutoffs for cells and genes
            min.cells = 0L,
            min.genes = 0L,
            # Default for UMI datasets
            is.expr = 0L,
            meta.data = as.data.frame(metrics(from))
        )

        # Check that the dimensions match exactly
        assert_are_identical(
            x = dim(from),
            y = dim(slot(to, "raw.data"))
        )

        # Stash metadata and rowRanges into `misc` slot
        bcbio <- list(
            "rowRanges" = rowRanges(from),
            "metadata" = metadata(from)
        )
        slot(to, "misc")[["bcbio"]] <- bcbio

        to
    }
)



#' @rdname coerce
#' @name coerce-seurat-SingleCellExperiment
#'
#' @section seurat to SingleCellExperiment:
#' Super basic S4 coercion support for taking the raw counts matrix from
#' a `seurat` class object and coercing to a `SingleCellExperiment`.
#'
#' @examples
#' # seurat to SingleCellExperiment ====
#' x <- as(seurat_small, "SingleCellExperiment")
#' class(x)
#' print(x)
setAs(
    from = "seurat",
    to = "SingleCellExperiment",
    function(from) {
        validObject(from)
        Convert(from = from, to = "sce")
    }
)
