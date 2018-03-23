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
#'
#' @return `data.frame` containing tSNE coordinates, sample metadata, and
#'   aggregate marker expression values (`mean`, `median`, and `sum`).
#'
#' @examples
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # seurat ====
#' genes <- counts(seurat_small) %>% rownames() %>% head()
#' print(genes)
#' fetchTSNEExpressionData(seurat_small, genes = genes) %>% glimpse()
NULL



# Constructors =================================================================
.fetchTSNEExpressionData.seurat <- function(  # nolint
    object,
    genes
) {
    tsne <- fetchTSNEData(object)

    data <- fetchGeneData(object, genes = genes)
    mean <- Matrix::rowMeans(data)
    median <- Biobase::rowMedians(data)
    sum <- Matrix::rowSums(data)

    cbind(tsne, mean, median, sum)
}



# Methods ======================================================================
#' @rdname fetchTSNEExpressionData
#' @export
setMethod(
    "fetchTSNEExpressionData",
    signature("seurat"),
    .fetchTSNEExpressionData.seurat
)
