## nocov start



#' @name defunct
#' @inherit AcidRoxygen::defunct description examples return seealso title
#' @inheritParams AcidRoxygen::params
#' @keywords internal
NULL



#' @name deprecated
#' @inherit AcidRoxygen::deprecated description examples return seealso title
#' @inheritParams AcidRoxygen::params
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



## v0.3.19 =====================================================================
#' @rdname deprecated
#' @export
plotGenesPerCell <- function(...) {
    .Deprecated("plotFeaturesPerCell")
    plotFeaturesPerCell(...)
}

#' @rdname deprecated
#' @export
plotUMIsPerCell <- function(...) {
    .Deprecated("plotCountsPerCell")
    plotCountsPerCell(...)
}

#' @rdname deprecated
#' @export
plotUMIsVsGenes <- function(...) {
    .Deprecated("plotCountsVsFeatures")
    plotCountsVsFeatures(...)
}



## nocov end
