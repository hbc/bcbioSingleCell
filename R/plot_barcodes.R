#' Plot all scRNA-Seq sample barcodes
#'
#' @author Michael Steinbaugh
#'
#' @param bcbio bcbio run object
#'
#' @export
plot_barcodes <- function(bcbio) {
    files <- list.files(bcbio$final_dir,
                        pattern = "*-barcodes.tsv",
                        recursive = TRUE,
                        include.dirs = TRUE,
                        full.names = TRUE)
    lapply(seq_along(files), function(a) {
        plot_barcode(files[a])
    }) %>% invisible
}
