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
#' fetchTSNEExpressionData(seurat, genes = genes) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname fetchTSNEExpressionData
#' @importFrom dplyr everything group_by select
#' @importFrom Matrix colMeans
#' @importFrom rlang !!! !! sym syms
#' @importFrom Seurat FetchData
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom tidyr gather
#' @export
setMethod(
    "fetchTSNEExpressionData",
    signature("seurat"),
    function(object, genes) {
        tsne <- fetchTSNEData(object)
        data <- FetchData(object, vars.all = genes) %>%
            as.data.frame()
        data[["geomean"]] <- Matrix::colMeans(t(data))
        data <- cbind(tsne, data)
        priorityCols <- c("gene", "cellID", "expression", "geomean")
        data %>%
            rownames_to_column("cellID") %>%
            gather(key = "gene",
                   value = "expression",
                   !!genes) %>%
            group_by(!!sym("gene")) %>%
            select(!!!syms(priorityCols), everything()) %>%
            arrange(!!!syms(priorityCols))
    })
