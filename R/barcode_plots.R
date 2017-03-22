#' Plot all scRNA-Seq sample barcodes
#'
#' @author Michael Steinbaugh
#'
#' @importFrom basejump bcbioProject
#'
#' @param `bcbio-nextgen` project list
#'
#' @export
barcode_plots <- function(project = NULL) {
    if (is.null(project)) {
        project <- basejump::bcbioProject("indrop_rnaseq")
    }
    files <- list.files(project$finalDir,
                        pattern = "*-barcodes.tsv",
                        recursive = TRUE,
                        include.dirs = TRUE,
                        full.names = TRUE)
    lapply(seq_along(files), function(a) {
        barcode_plot(files[a])
    }) %>% invisible
}
