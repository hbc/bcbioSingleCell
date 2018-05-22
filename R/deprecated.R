# nocov start



#' Deprecated Functions
#'
#' @name deprecated
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return [.Defunct()] or [.Deprecated()].
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
#' @export
darkTheme <- function(...) {
    .Deprecated("basejump::theme_midnight")
    basejump::theme_midnight(...)
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
    .Deprecated("bcbioBase::plotQuantileHeatmap")
    bcbioBase::plotQuantileHeatmap(...)
}



# v0.0.25 ======================================================================
#' @rdname deprecated
#' @export
plotKnownMarkers <- function(object, knownMarkers, ...) {
    .Deprecated("plotKnownMarkersDetected")
    plotKnownMarkersDetected(
        object,
        knownMarkersDetected = knownMarkers,
        ...
    )
}



#' @rdname deprecated
#' @export
readMarkers <- function(...) {
    .Deprecated("readCellTypeMarkers")
    readCellTypeMarkers(...)
}



# v0.0.26 ======================================================================
#' @rdname deprecated
#' @export
readMarkersFile <- function(...) {
    .Deprecated("readCellTypeMarkers")
    readCellTypeMarkers(...)
}



# v0.1.0 =======================================================================
#' @rdname deprecated
#' @export
calculateMetrics <- function(...) {
    .Deprecated("metrics")
    metrics(...)
}



#' @rdname deprecated
#' @export
readCellTypeMarkersFile <- function(...) {
    .Deprecated("readCellTypeMarkers")
    readCellTypeMarkers(...)
}



# v0.1.1 =======================================================================
#' @rdname deprecated
#' @export
inflectionPoint <- function(...) {
    .Defunct("DropletUtils::barcodeRanks")
}



#' @rdname deprecated
#' @export
plotCumulativeUMIsPerCell <- function(...) {
    .Defunct("plotUMIsPerCell")
}



# v0.1.2 =======================================================================
#' @rdname deprecated
#' @export
loadCellRanger <- function(...) {
    .Deprecated("readCellRanger")
    readCellRanger(...)
}



#' @rdname deprecated
#' @export
loadSingleCell <- function(...) {
    .Deprecated("bcbioSingleCell")
    bcbioSingleCell(...)
}



# v0.1.10 ======================================================================
#' @rdname deprecated
#' @export
plotMarker <- function(object, gene, ...) {
    .Deprecated("plotMarkers")
    plotMarkers(object, genes = gene, ...)
}



# nocov end
