#' Prepare Single-Cell RNA-Seq RMarkdown Template
#'
#' @rdname prepareSingleCellTemplate
#' @name prepareSingleCellTemplate
#' @family RMarkdown Utilities
#' @author Michael Steinbaugh
#'
#' @return No value.
#'
#' @examples
#' \dontrun{
#' prepareSingleCellTemplate()
#' }
NULL



# Methods ====
#' @rdname prepareSingleCellTemplate
#' @export
setMethod("prepareSingleCellTemplate", "missing", function(object) {
    prepareTemplate(
        system.file("rmarkdown/shared", package = "bcbioSingleCell"))
})
