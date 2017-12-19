#' Fetch t-SNE Locations and Cellular Metadata
#'
#' @rdname fetchTSNEData
#' @name fetchTSNEData
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return [data.frame] of t-SNE points and metadata for each cell.
#'
#' @seealso [Seurat::TSNEPlot()].
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # seurat
#' fetchTSNEData(seurat) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname fetchTSNEData
#' @export
setMethod(
    "fetchTSNEData",
    signature("seurat"),
    function(object) {
        .fetchDimDataSeurat(
            object,
            dimCode = c(x = "tSNE_1", y = "tSNE_2")
        )
    })
