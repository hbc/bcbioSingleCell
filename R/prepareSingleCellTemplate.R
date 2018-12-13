#' Prepare Single-Cell RNA-Seq R Markdown Template
#'
#' @name prepareSingleCellTemplate
#' @author Michael Steinbaugh
#' @inherit basejump::prepareTemplate
#' @export
#'
#' @examples
#' \dontrun{
#' ## prepareSingleCellTemplate()
#' }
prepareSingleCellTemplate <- function(overwrite = FALSE) {
    prepareTemplate(package = "bcbioSingleCell", overwrite = overwrite)
}
