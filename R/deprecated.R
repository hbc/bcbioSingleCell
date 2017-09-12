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
plotClusters <- function(...) {
    .Deprecated("plotMarkers")
    plotMarkers(...)
}
