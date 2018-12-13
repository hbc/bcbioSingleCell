#' Prepare Single-Cell RNA-Seq R Markdown Template
#'
#' @name prepareSingleCellTemplate
#' @author Michael Steinbaugh
#' @inherit basejump::prepareTemplate
#' @export
#'
#' @examples
#' x <- prepareSingleCellTemplate()
#' x
#'
#' ## Clean up.
#' unlink(c(
#'     "_footer.Rmd",
#'     "_header.Rmd",
#'     "_output.yaml",
#'     "_setup.R",
#'     "bibliography.bib"
#' ))
prepareSingleCellTemplate <- function(overwrite = FALSE) {
    prepareTemplate(package = "bcbioSingleCell", overwrite = overwrite)
}
