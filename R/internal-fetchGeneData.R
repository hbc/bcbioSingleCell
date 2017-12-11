#' Fetch Gene Expression Data from Seurat
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom Seurat FetchData
#'
#' @param object seurat.
#' @param genes Gene identifiers (matrix rownames).
.fetchGeneDataSeurat <- function(object, genes) {
    data <- Seurat::FetchData(object, vars.all = genes)
    data <- as.data.frame(data)
    if (!identical(as.character(genes), colnames(data))) {
        stop("'Seurat::FetchData()' return doesn't match 'genes' input")
    }
    data
}
