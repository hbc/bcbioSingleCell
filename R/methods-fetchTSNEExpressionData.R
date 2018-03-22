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
#' @inheritParams general
#'
#' @param genes Genes identifiers (matching the rownames in the object),
#'   of which to get expression data.
#'
#' @return [data.frame] containing tSNE coordinates, sample metadata, and
#'   aggregate marker expression values (`mean`, `median`, and `sum`).
#'
#' @examples
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' genes <- counts(seurat) %>% rownames() %>% head()
#' print(genes)
#'
#' fetchTSNEExpressionData(seurat, genes = genes) %>%
#'     glimpse()
NULL



# Constructors =================================================================
.fetchTSNEExpressionData.seurat <- function(  # nolint
    object,
    genes) {
    tsne <- fetchTSNEData(object)

    # Gene aggregate math
    data <- fetchGeneData(object, genes = genes)
    # TODO Okay to just use S4Vectors methods here?
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
    .fetchTSNEExpressionData.seurat)
