#' Prepare Single-Cell RNA-Seq RMarkdown Template
#'
#' @name prepareSingleCellTemplate
#' @family R Markdown Functions
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase prepareTemplate
#'
#' @inheritParams general
#'
#' @return No value.
#' @export
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
prepareSingleCellTemplate <- function() {
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
