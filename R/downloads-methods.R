#' File Downloads
#'
#' If file isn't present, download latest version from the
#' [HBC website](http://bioinformatics.sph.harvard.edu/).
#'
#' File download utility functions for RMarkdown knit reports.
#'
#' @rdname downloads
#' @name downloads
#' @author Michael Steinbaugh
#'
#' @return No return.
#'
#' @examples
#' downloads()
#' downloads("setup.R")
NULL



# Constructors ====
.downloads <- function(file) {
    sapply(seq_along(file), function(a) {
        if (!file.exists(file[[a]])) {
            download.file(
                file.path("http://bioinformatics.sph.harvard.edu",
                          "bcbioSinglecell", "downloads", file[[a]]),
                destfile = file[[a]])
        }
    }) %>% invisible
}



# Methods ====
#' @rdname downloads
#' @export
setMethod("downloads", "missing", function(object) {
    .downloads(
        c("_output.yaml",
          "_footer.Rmd",
          "_header.Rmd",
          "bcbioSinglecell.bib",
          "load_run.R",
          "setup.R"))
})



#' @rdname downloads
#' @export
setMethod("downloads", "character", function(object) {
    .downloads(object)
})
