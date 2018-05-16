#' Prepare Single-Cell RNA-Seq R Markdown Template
#'
#' @name prepareSingleCellTemplate
#' @family R Markdown Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase prepareTemplate
#'
#' @inherit bcbioBase::prepareTemplate
#'
#' @export
#'
#' @examples
#' x <- prepareSingleCellTemplate()
#' x
#'
#' # Clean up
#' unlink(c(
#'     "_footer.Rmd",
#'     "_header.Rmd",
#'     "_output.yaml",
#'     "_setup.R",
#'     "bibliography.bib"
#' ))
prepareSingleCellTemplate <- function(overwrite = FALSE) {
    prepareTemplate(
        file = c(
            "_footer.Rmd",
            "_header.Rmd",
            "_output.yaml",
            "_setup.R",
            "bibliography.bib"
        ),
        sourceDir = system.file(
            "rmarkdown/shared",
            package = "bcbioSingleCell"
        ),
        overwrite = overwrite
    )
}
