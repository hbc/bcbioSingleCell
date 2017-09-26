#' Fetch t-SNE Locations, Cellular Metadata and Expression of Genes
#'
#' This gets t-SNE locations, cellular metadata, expression of genes and
#' the geometric mean of the gene expression from a [seurat] object.
#'
#' @rdname fetchTSNEExpressionData
#' @name fetchTSNEExpressionData
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param genes Genes (by symbol name) of which to get expression data.
#'
#' @return tidy [data.frame] of t-SNE points, cellular metadata, and gene
#'   expression.
#'
#' @examples
#' \dontrun{
#' data(seurat)
#' genes <- head(rownames(seurat@raw.data))
#' tsne <- fetchTSNEExpressionData(seurat, genes)
#' }
NULL



# Methods ====
#' @rdname fetchTSNEExpressionData
#' @export
setMethod("fetchTSNEExpressionData", "seurat", function(
    object, genes) {
    dat <- FetchData(object, vars.all = genes) %>%
        as.data.frame()
    dat[["geomean"]] <- colMeans(t(dat))
    dat <- dat %>%
        rownames_to_column("cell")
    fetchTSNEData(object) %>%
        left_join(dat, by = "cell") %>%
        gather(key = "gene",
               value = "expression",
               !!genes)
})
