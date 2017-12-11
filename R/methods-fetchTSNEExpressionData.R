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
#' @param format Gene identifier format. Supports `ensgene` or `symbol`.
#'
#' @return [tibble] grouped by `gene`, containing t-SNE points, cellular
#'   metadata, and gene expression.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' symbol <- counts(seurat) %>% rownames() %>% head()
#' print(symbol)
#'
#' ensgene <- bcbio(seurat, "gene2symbol") %>%
#'     .[which(.[["symbol"]] %in% symbol), "ensgene", drop = TRUE]
#' print(ensgene)
#'
#' fetchTSNEExpressionData(
#'     seurat,
#'     genes = symbol,
#'     format = "symbol") %>%
#'     glimpse()
#' fetchTSNEExpressionData(
#'     seurat,
#'     genes = ensgene,
#'     format = "ensgene") %>%
#'     glimpse()
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
    function(
        object,
        genes,
        format = "symbol") {
        priorityCols <- c("gene", "cellID", "expression", "geomean")
        .checkFormat(format)
        if (format == "ensgene") {
            genes <- .convertGenesToSymbols(object, genes = genes)
        }
        data <- .fetchGeneDataSeurat(object, genes = genes)
        data[["geomean"]] <- Matrix::colMeans(t(data))
        tsne <- fetchTSNEData(object)
        cbind(tsne, data) %>%
            rownames_to_column("cellID") %>%
            gather(key = "gene",
                   value = "expression",
                   !!genes) %>%
            group_by(!!sym("gene")) %>%
            select(!!!syms(priorityCols), everything()) %>%
            arrange(!!!syms(priorityCols))
    })
