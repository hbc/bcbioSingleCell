#' Fetch Gene Expression Data from Seurat
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom Seurat FetchData
#'
#' @param object [seurat].
#' @param genes Gene identifiers (matrix rownames).
#'
#' @return [matrix].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' genes <- slot(seurat, "data") %>%
#'     rownames() %>%
#'     head()
#'
#' .fetchGeneData.seurat(seurat, genes = genes) %>%
#'     glimpse()
.fetchGeneData.seurat <- function(object, genes) {  # nolint
    data <- Seurat::FetchData(object, vars.all = genes)
    if (!identical(
        as.character(genes),
        as.character(colnames(data))
    )) {
        abort("`Seurat::FetchData()` return doesn't match `genes` input")
    }
    data
}
