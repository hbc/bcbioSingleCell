#' File Downloads
#'
#' If file isn't present, download latest version from the
#' [HBC website](http://bioinformatics.sph.harvard.edu/).
#'
#' File download utility functions for RMarkdown knit reports.
#'
#' @rdname download
#' @name download
#'
#' @return No return.
#'
#' @examples
#' download()
#' download("setup.R")
NULL



# Constructors ====
.download <- function(file) {
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
#' @rdname download
#' @export
setMethod("download", "missing", function() {
    .download(
        c("_output.yaml",
          "_footer.Rmd",
          "_header.Rmd",
          "bcbioSinglecell.bib",
          "setup.R"))
})



#' @rdname download
#' @export
setMethod("download", "character", function(object) {
    .download(object)
})
