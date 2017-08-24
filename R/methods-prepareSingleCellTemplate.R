#' Prepare Single-Cell RNA-Seq RMarkdown Template
#'
#' @rdname prepareSingleCellTemplate
#' @name prepareSingleCellTemplate
#'
#' @return No value.
NULL



# Methods ====
#' @rdname prepareSingleCellTemplate
#' @export
setMethod("prepareSingleCellTemplate", "missing", function(object) {
    prepareTemplate(package = "bcbioSinglecell")
})
