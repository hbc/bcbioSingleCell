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



# v0.1.0 =======================================================================
#' @rdname deprecated
#' @export
calculateMetrics <- function(...) {
    .Deprecated("metrics")
    metrics(...)
}



# v0.1.1 =======================================================================
#' @rdname deprecated
#' @export
inflectionPoint <- function(...) {
    .Defunct("plotBarcodeRanks")
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



# nocov end
