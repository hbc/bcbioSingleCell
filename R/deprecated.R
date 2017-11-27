#' Deprecated Functions
#'
#' @rdname deprecated
#' @name deprecated
#' @keywords internal
#'
#' @inheritParams AllGenerics
#'
#' @return Soft deprecation to new functions.
NULL



# 0.0.18 ====
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



# 0.0.19 ====
#' @rdname deprecated
#' @export
loadSingleCellRun <- function(...) {
    .Deprecated("loadSingleCell")
    loadSingleCell(...)
}



# 0.0.23 ====
#' @rdname deprecated
#' @export
plotFeatures <- function(object, features, ...) {
    .Deprecated("plotFeatureTSNE")
    plotFeatureTSNE(object = object, feature = features, ...)
}
