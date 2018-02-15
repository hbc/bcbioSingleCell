#' Prepare Single-Cell RNA-Seq RMarkdown Template
#'
#' @rdname prepareSingleCellTemplate
#' @name prepareSingleCellTemplate
#' @family RMarkdown Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return No value.
#'
#' @examples
#' files <- c(
#'     "_footer.Rmd",
#'     "_header.Rmd",
#'     "_output.yaml",
#'     "_setup.R",
#'     "bibliography.bib"
#' )
#' prepareSingleCellTemplate()
#' all(file.exists(files))
#' unlink(files)
NULL



# Constructors =================================================================
.prepareSingleCellTemplate <- function(object) {
    prepareTemplate(
        c(
            "_footer.Rmd",
            "_header.Rmd",
            "_output.yaml",
            "_setup.R",
            "bibliography.bib"
        ),
        sourceDir = system.file(
            "rmarkdown/shared",
            package = "bcbioSingleCell"
        )
    )
}



# Methods ======================================================================
#' @rdname prepareSingleCellTemplate
#' @importFrom bcbioBase prepareTemplate
#' @export
setMethod(
    "prepareSingleCellTemplate",
    signature("missing"),
    .prepareSingleCellTemplate)
