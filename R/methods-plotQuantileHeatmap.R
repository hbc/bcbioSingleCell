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
#' @inheritParams general
#'
#' @param n The number of breaks to create.
#' @param annotationCol *Optional.* [data.frame] that defines annotation
#'   mappings for the columns.
#' @param clusterCols Logical determining if columns should be arranged with
#'   hierarchical clustering. Alternatively, can define an `hclust` object.
#' @param clusterRows Logical determining if rows should be arranged with
#'   hierarchical clustering. Alternatively, can define an `hclust` object.
#' @param color Colors to use for plot. Defaults to the [viridis::viridis()]
#'   palette.
#' @param legendColor Colors to use for legend labels. Defaults to the
#'   [viridis::viridis()] palette.
#' @param title *Optional.* Plot title.
#'
#' @examples
#' load(system.file("extdata/filtered.rda", package = "bcbioSingleCell"))
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
        n = 5L,
        annotationCol = NA,
        clusterCols = FALSE,
        clusterRows = FALSE,
        color = viridis::viridis,
        legendColor = viridis::viridis,
        title = NULL) {
        counts <- counts(object)
        plotQuantileHeatmap(
            object = counts,
            n = n,
            annotationCol = annotationCol,
            clusterCols = clusterCols,
            clusterRows = clusterRows,
            color = color,
            legendColor = legendColor,
            title = title)
    })



#' @rdname plotQuantileHeatmap
#' @importFrom viridis viridis
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("seurat"),
    function(
        object,
        n = 5L,
        annotationCol = NA,
        clusterCols = FALSE,
        clusterRows = FALSE,
        color = viridis::viridis,
        legendColor = viridis::viridis,
        title = NULL) {
        counts <- counts(object, normalized = FALSE)
        plotQuantileHeatmap(
            object = counts,
            n = n,
            annotationCol = annotationCol,
            clusterCols = clusterCols,
            clusterRows = clusterRows,
            color = color,
            legendColor = legendColor,
            title = title)
    })
