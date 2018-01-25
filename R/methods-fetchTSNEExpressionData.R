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
#' @param genes Genes identifiers (matching the rownames in the object),
#'   of which to get expression data.
#'
#' @return [tibble] grouped by `gene`, containing t-SNE points, cellular
#'   metadata, and gene expression.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' genes <- counts(seurat) %>% rownames() %>% head()
#' print(genes)
#'
#' fetchTSNEExpressionData(seurat, genes = genes) %>%
#'     glimpse()
NULL



# Constructors =================================================================
#' @importFrom dplyr everything group_by select
#' @importFrom rlang !!! !! sym syms
#' @importFrom Seurat FetchData
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tidyr gather
.fetchTSNEExpressionData.seurat <- function(  # nolint
    object,
    genes) {
    tsne <- fetchTSNEData(object)

    # Gene aggregate math
    data <- fetchGeneData(object, genes = genes)
    mean <- rowMeans(data)
    median <- rowMedians(data)
    sum <- rowSums(data)

    cbind(tsne, mean, median, sum)
}



# Methods ======================================================================
#' @rdname fetchTSNEExpressionData
#' @export
setMethod(
    "fetchTSNEExpressionData",
    signature("seurat"),
    .fetchTSNEExpressionData.seurat)
