#' Deprecated Functions
#'
#' @rdname deprecated
#' @name deprecated
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return Soft deprecation to new functions.
NULL



# v0.0.18 ======================================================================
#' @rdname deprecated
#' @export
plotClusters <- function(...) {
    .Deprecated("plotMarkers")
    plotMarkers(...)
}



#' @rdname deprecated
#' @export
plotTSNEExpressionData <- function(...) {
    .Deprecated("plotMarkerTSNE")
    plotMarkerTSNE(...)
}



# v0.0.19 ======================================================================
#' @rdname deprecated
#' @export
loadSingleCellRun <- function(...) {
    .Deprecated("loadSingleCell")
    loadSingleCell(...)
}



# v0.0.23 ======================================================================
#' @rdname deprecated
#' @export
plotFeatures <- function(object, features, ...) {
    .Deprecated("plotFeatureTSNE")
    plotFeatureTSNE(object = object, feature = features, ...)
}



# v0.0.24 ======================================================================
#' @rdname deprecated
#' @importFrom basejump midnightTheme
#' @export
darkTheme <- function(...) {
    .Deprecated("midnightTheme")
    midnightTheme(...)
}



#' @rdname deprecated
#' @export
pcCutoff <- function(...) {
    .Deprecated("plotPCElbow")
    plotPCElbow(...)
}



#' @rdname deprecated
#' @export
quantileHeatmap <- function(...) {
    .Deprecated("plotQuantileHeatmap")
    plotQuantileHeatmap(...)
}



# v0.0.25 ======================================================================
#' @rdname deprecated
#' @export
plotKnownMarkers <- function(
    object,
    knownMarkers,
    ...) {
    .Deprecated("plotKnownMarkersDetected")
    plotKnownMarkersDetected(
        object,
        knownMarkersDetected = knownMarkers,
        ...)
}



#' @rdname deprecated
#' @export
readMarkers <- function(...) {
    .Deprecated("readCellTypeMarkersFile")
    readCellTypeMarkersFile(...)
}



# v0.0.26 ======================================================================
readMarkersFile <- function(...) {
    .Deprecated("readCellTypeMarkersFile")
    readCellTypeMarkersFile(...)
}
