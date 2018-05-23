# nocov start



#' Deprecated Functions
#'
#' @name deprecated
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return [.Deprecated()].
NULL



#' Defunct Functions
#'
#' @name defunct
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return [.Defunct()].
NULL



# v0.0.18 ======================================================================
#' @rdname defunct
#' @export
plotClusters <- function(...) {
    .Defunct("plotMarkers")
}



#' @rdname defunct
#' @export
plotTSNEExpressionData <- function(...) {
    .Defunct("plotMarkerTSNE")
}



# v0.0.19 ======================================================================
#' @rdname defunct
#' @export
loadSingleCellRun <- function(...) {
    .Defunct("bcbioSingleCell")
}



# v0.0.24 ======================================================================
#' @rdname defunct
#' @export
darkTheme <- function(...) {
    .Defunct("basejump::theme_midnight")
}



#' @rdname defunct
#' @export
pcCutoff <- function(...) {
    .Defunct("plotPCElbow")
}



#' @rdname defunct
#' @export
quantileHeatmap <- function(...) {
    .Defunct("plotQuantileHeatmap")
}



# v0.0.25 ======================================================================
#' @rdname defunct
#' @export
plotKnownMarkers <- function(...) {
    .Defunct("plotKnownMarkersDetected")
}



#' @rdname defunct
#' @export
readMarkers <- function(...) {
    .Defunct("readCellTypeMarkers")
}



# v0.0.26 ======================================================================
#' @rdname defunct
#' @export
readMarkersFile <- function(...) {
    .Defunct("readCellTypeMarkers")
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
#' @rdname defunct
#' @export
plotMarker <- function(...) {
    .Defunct("plotMarkers")
}



# nocov end
