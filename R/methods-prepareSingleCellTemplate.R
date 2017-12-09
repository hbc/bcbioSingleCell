#' Prepare Single-Cell RNA-Seq RMarkdown Template
#'
#' @rdname prepareSingleCellTemplate
#' @name prepareSingleCellTemplate
#' @family RMarkdown Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return No value.
#'
#' @examples
#' \dontrun{
#' prepareSingleCellTemplate()
#' }
NULL



# Methods ======================================================================
#' @rdname prepareSingleCellTemplate
#' @export
setMethod(
    "prepareSingleCellTemplate",
    signature("missing"),
    function(object) {
        prepareTemplate(
            sourceDir = system.file(
                "rmarkdown/shared",
                package = "bcbioSingleCell"))
    })
