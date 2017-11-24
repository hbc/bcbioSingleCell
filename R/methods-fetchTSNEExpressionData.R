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
#'
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
#' @importFrom dplyr left_join
#' @importFrom Matrix colMeans
#' @importFrom Seurat FetchData
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @export
setMethod(
    "fetchTSNEExpressionData",
    signature("seurat"),
    function(object, genes) {
        data <- FetchData(object, vars.all = genes)
        data <- as.data.frame(data)
        data[["geomean"]] <- Matrix::colMeans(t(data))
        data <- rownames_to_column(data, "cell")
        fetchTSNEData(object) %>%
            left_join(data, by = "cell") %>%
            gather(key = "gene",
                   value = "expression",
                   !!genes)
    })
