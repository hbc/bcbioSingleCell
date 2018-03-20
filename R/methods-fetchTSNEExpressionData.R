#' Fetch t-SNE Locations, Cellular Metadata and Expression of Genes
#'
#' This gets t-SNE locations, cellular metadata, expression of genes and
#' the geometric mean of the gene expression from a [seurat] object.
#'
#' @name fetchTSNEExpressionData
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#' @param genes Genes identifiers (matching the rownames in the object),
#'   of which to get expression data.
#'
#' @return `data.frame` containing tSNE coordinates, sample metadata, and
#'   aggregate marker expression values (`mean`, `median`, and `sum`).
#'
#' @examples
#' # seurat ====
#' genes <- head(rownames(pbmc_small))
#' print(genes)
#' fetchTSNEExpressionData(pbmc_small, genes = genes) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname fetchTSNEExpressionData
#' @importFrom Biobase rowMedians
#' @export
setMethod(
    "fetchTSNEExpressionData",
    signature("seurat"),
    function(object, genes) {
        assert_is_character(genes)
        tsne <- fetchTSNEData(object)

        # Gene aggregate math
        data <- fetchGeneData(object, genes = genes)
        mean <- Matrix::rowMeans(data)
        median <- rowMedians(data)
        sum <- Matrix::rowSums(data)

        cbind(tsne, mean, median, sum)
    }
)
