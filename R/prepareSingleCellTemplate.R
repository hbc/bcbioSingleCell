#' Prepare Single-Cell RNA-Seq R Markdown Template
#'
#' @name prepareSingleCellTemplate
#' @family R Markdown Functions
#' @author Michael Steinbaugh
#'
#' @inherit basejump::prepareTemplate
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
    package <- "bcbioSingleCell"
    prepareTemplate(
        package = package,
        sourceDir = system.file(
            "rmarkdown/shared",
            package = package,
            mustWork = TRUE
        ),
        overwrite = overwrite
    )
}
