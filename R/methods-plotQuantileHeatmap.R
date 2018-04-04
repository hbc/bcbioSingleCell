#' Plot Heatmap with Quantile Breaks
#'
#' @details This is helpful for visualizing single cell data.
#'   Ideas and code from: http://slowkow.com/notes/heatmap-tutorial/
#'
#' @name plotQuantileHeatmap
#' @family Quality Control Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @importFrom bcbioBase plotQuantileHeatmap
#'
#' @inherit bcbioBase::plotQuantileHeatmap
#'
#' @examples
#' # bcbioSingleCell ====
#' plotQuantileHeatmap(
#'     object = bcb_small,
#'     showColnames = FALSE,
#'     showRownames = FALSE
#' )
#'
#' # seurat ====
#' plotQuantileHeatmap(
#'     object = seurat_small,
#'     showColnames = FALSE,
#'     showRownames = FALSE
#' )
NULL



# Methods ======================================================================
#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("seurat"),
    getMethod("plotQuantileHeatmap", "SummarizedExperiment")
)
