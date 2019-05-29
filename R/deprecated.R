# nocov start



#' @name defunct
#' @inherit basejump::defunct
#' @keywords internal
NULL

#' @name deprecated
#' @inherit basejump::deprecated
#' @keywords internal
NULL



# v0.2.2 =======================================================================
#' @rdname defunct
#' @export
readCellRanger <- function(...) {
    .Defunct(msg = "CellRanger support has migrated to Chromium package.")
}



# v0.3.12 ======================================================================
#' @rdname deprecated
#' @export
prepareSingleCellTemplate <- function(...) {
    .Deprecated("prepareTemplate(package = \"bcbioSingleCell\")")
    prepareTemplate(package = "bcbioSingleCell", ...)
}



# nocov end
