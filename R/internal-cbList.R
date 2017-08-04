#' Cellular Barcodes List
#'
#' @rdname internal-cbList
#' @keywords internal
#'
#' @param sampleDirs Sample directories.
#'
#' @return [list].
.cbList <- function(sampleDirs) {
    files <- sampleDirs %>%
        file.path(paste(basename(.), "barcodes.tsv", sep = "-")) %>%
        set_names(names(sampleDirs))
    if (!all(file.exists(files))) {
        stop("Cellular barcode file missing")
    }
    message("Reading cellular barcode distributions")
    pblapply(seq_along(files), function(a) {
        .readCBFile(files[a])
    }) %>% set_names(names(sampleDirs))
}
