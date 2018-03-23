#' Fetch Gene Expression Data
#'
#' @name fetchGeneData
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `matrix`.
#'
#' @examples
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # seurat
#' genes <- slot(seurat, "data") %>%
#'     rownames() %>%
#'     head()
#' fetchGeneData(seurat, genes = genes) %>%
#'     head()
NULL



# Constructors =================================================================
.fetchGeneData.seurat <- function(object, genes) {  # nolint
    data <- FetchData(object, vars.all = genes)
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
