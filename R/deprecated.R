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



#' @rdname deprecated
#' @export
metadataTable <- function(...) {
    .Deprecated("sampleMetadata")
    sampleMetadata(...)
}



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
