#' Fetch Gene Expression Data from Seurat
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom Seurat FetchData
#'
#' @param object [seurat].
#' @param genes Gene identifiers (matrix rownames).
.fetchGeneData.Seurat <- function(object, genes) {
    data <- Seurat::FetchData(object, vars.all = genes)
    if (!identical(
        as.character(genes),
        as.character(colnames(data))
    )) {
        stop("'Seurat::FetchData()' return doesn't match 'genes' input")
    }
    data
}
