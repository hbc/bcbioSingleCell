#' Plot Heatmap with Quantile Breaks
#'
#' @details This is helpful for more usefully visualizing single cell data.
#'   Ideas and code from: http://slowkow.com/notes/heatmap-tutorial/
#'
#' @rdname plotQuantileHeatmap
#' @name plotQuantileHeatmap
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @importFrom basejump plotQuantileHeatmap
#'
#' @inherit basejump::plotQuantileHeatmap
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "filtered.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' plotQuantileHeatmap(filtered)
NULL



# Methods ======================================================================
#' @rdname plotQuantileHeatmap
#' @importFrom viridis viridis
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("bcbioSingleCell"),
    function(
        object,
        n = 10,
        annotationCol = NA,
        clusterRows = FALSE,
        clusterCols = FALSE,
        color = viridis::viridis) {
        counts <- counts(object)
        plotQuantileHeatmap(
            object = counts,
            n = n,
            annotationCol = annotationCol,
            clusterRows = clusterRows,
            clusterCols = clusterCols,
            color = color)
    })



#' @rdname plotQuantileHeatmap
#' @importFrom viridis viridis
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("seurat"),
    function(
        object,
        n = 10,
        annotationCol = NA,
        clusterRows = FALSE,
        clusterCols = FALSE,
        color = viridis::viridis) {
        counts <- counts(object, normalized = FALSE)
        plotQuantileHeatmap(
            object = counts,
            n = n,
            annotationCol = annotationCol,
            clusterRows = clusterRows,
            clusterCols = clusterCols,
            color = color)
    })
