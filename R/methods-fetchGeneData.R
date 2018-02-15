#' Fetch Gene Expression Data
#'
#' @rdname fetchGeneData
#' @name fetchGeneData
#'
#' @inheritParams general
#'
#' @param genes Gene identifiers (matrix rownames).
#'
#' @return [matrix].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # seurat
#' genes <- slot(seurat, "data") %>%
#'     rownames() %>%
#'     head()
#' fetchGeneData(seurat, genes = genes) %>%
#'     head()
NULL



# Constructors =================================================================
#' @importFrom Seurat FetchData
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



# Methods ======================================================================
#' @rdname fetchGeneData
#' @export
setMethod(
    "fetchGeneData",
    signature("seurat"),
    .fetchGeneData.seurat)
