#' Fetch t-SNE locations, cellular metadata and expression of genes
#'
#' This gets t-SNE locations, cellular metadata, expression of genes and
#' the geometric mean of the gene expression from a seurat object.
#'
#' @rdname fetchTSNEExpressionData
#' @name fetchTSNEExpressionData
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @param genes Genes (by symbol name) of which to get expression data.
#' @return tidy [data.frame] of t-SNE points, cellular metadata, and gene
#'   expression.



# Methods ====
#' @rdname fetchTSNEExpressionData
#' @export
setMethod("fetchTSNEExpressionData", "seurat", function(
    object, genes) {
    dat <- FetchData(object, vars.all = genes) %>%
        as.data.frame
    dat[["geomean"]] <- colMeans(t(dat))
    dat <- dat %>%
        rownames_to_column("cell")
    fetchTSNEData(object) %>%
        left_join(dat, by = "cell") %>%
        gather(gene, expression, !!genes)
})
