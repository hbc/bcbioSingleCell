#' Methods for Coercing an Object to a Class
#'
#' @name coerce
#' @aliases as
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
#' @name coerce-bcbioSingleCell-seurat
#' @section bcbioSingleCell to seurat:
#' Interally [Seurat::CreateSeuratObject()] is called without applying any
#' additional filtering cutoffs, since we have already defined them during our
#' quality control analysis. Here we are passing the raw gene-level counts of
#' the filtered cells into a new `seurat` class object.
setAs(
    from = "bcbioSingleCell",
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

        # Require filtering of cells and genes
        from <- .applyFilterCutoffs(from)

        # Prepare counts matrix with gene symbols as rownames
        counts <- counts(from)
        g2s <- gene2symbol(from)
        assertIsGene2symbol(g2s)
        rownames <- make.names(g2s[["geneName"]], unique = TRUE)
        stopifnot(!any(duplicated(rownames)))
        names(rownames) <- g2s[["geneID"]]
        rownames(counts) <- rownames

        # New seurat object
        seurat <- Seurat::CreateSeuratObject(
            raw.data = counts,
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
            y = dim(slot(seurat, "raw.data"))
        )

        # Stash metadata and rowRanges into `misc` slot
        bcbio <- list(
            "rowRanges" = rowRanges(from),
            "metadata" = metadata(from)
        )
        slot(seurat, "misc")[["bcbio"]] <- bcbio

        seurat
    }
)
