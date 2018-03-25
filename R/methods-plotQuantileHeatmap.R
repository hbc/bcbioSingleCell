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
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' plotQuantileHeatmap(bcb_small)
#'
#' # seurat ====
#' plotQuantileHeatmap(seurat_small)
NULL



# Methods ======================================================================
#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("seurat"),
    getMethod("plotQuantileHeatmap", "SummarizedExperiment")
)
