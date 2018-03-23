#' Fetch t-SNE Locations and Cellular Metadata
#'
#' @name fetchTSNEData
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `data.frame` of t-SNE points and metadata for each cell.
#'
#' @seealso [Seurat::TSNEPlot()].
#'
#' @examples
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
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
        .fetchDRData.seurat(
            object,
            dimCode = c(x = "tSNE_1", y = "tSNE_2")
        )
    })
