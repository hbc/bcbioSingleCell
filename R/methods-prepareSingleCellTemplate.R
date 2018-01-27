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
#' files <- c(
#'     "_footer.Rmd", "_header.Rmd", "_output.yaml",
#'     "bibliography.bib", "setup.R")
#' prepareSingleCellTemplate()
#' all(file.exists(files))
#' unlink(files)
NULL



# Methods ======================================================================
#' @rdname prepareSingleCellTemplate
#' @importFrom bcbioBase prepareTemplate
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
