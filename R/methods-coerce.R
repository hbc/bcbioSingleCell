#' Methods for Coercing an Object to a Class
#'
#' @name coerce
#' @aliases as
#' @family S4 Class Definition
#' @author Michael Steinbaugh
#'
#' @inherit bcbioBase::coerce
#'
#' @seealso
#' - [methods::as()].
#' - [methods::coerce()].
#'
#' @examples
#' # SingleCellExperiment to seurat ====
#' x <- as(bcb_small, "seurat")
#' class(x)
#' print(x)
NULL



# Methods ======================================================================
#' @rdname coerce
#' @name coerce-SingleCellExperiment-seurat
#' @section SingleCellExperiment to seurat:
#' Interally [Seurat::CreateSeuratObject()] is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new `seurat` class object.
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
        rownames <- rownames(from)
        from <- convertGenesToSymbols(from)

        # Ensure that genes are unique valid names.
        # Note that any "-" in gene names will be sanitized to "." here.
        rownames(from) <- make.names(rownames(from), unique = TRUE)

        # Create the seurat object
        to <- Seurat::CreateSeuratObject(
            raw.data = counts(from),
            project = "bcbioSingleCell",
            # Already applied filtering cutoffs for cells and genes
            min.cells = 0L,
            min.genes = 0L,
            # Default for UMI datasets
            is.expr = 0L,
            meta.data = metrics(from)
        )

        # Check that the dimensions match exactly
        assert_are_identical(
            x = dim(from),
            y = dim(slot(to, "raw.data"))
        )

        # Stash metadata and rowRanges into `misc` slot
        bcbio <- list(
            "rownames" = rownames,
            "rowRanges" = rowRanges(from),
            "metadata" = metadata(from)
        )
        slot(to, "misc")[["bcbio"]] <- bcbio

        to
    }
)
