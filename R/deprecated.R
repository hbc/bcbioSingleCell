## nocov start



#' @name defunct
#' @inherit acidroxygen::defunct description examples return seealso title
#' @inheritParams acidroxygen::params
#' @keywords internal
NULL



#' @name deprecated
#' @inherit acidroxygen::deprecated description examples return seealso title
#' @inheritParams acidroxygen::params
#' @keywords internal
NULL



## v0.2.2 ======================================================================
#' @rdname defunct
#' @export
readCellRanger <- function(...) {
    .Defunct(msg = "CellRanger support has migrated to Chromium package.")
}



## v0.3.12 =====================================================================
#' @rdname defunct
#' @export
prepareSingleCellTemplate <- function(...) {
    .Defunct("prepareTemplate(package = \"bcbioSingleCell\")")
}



## nocov end
