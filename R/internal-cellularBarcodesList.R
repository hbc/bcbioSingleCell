#' Cellular Barcodes List
#'
#' @author Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom pbapply pblapply
#' @importFrom stats setNames
#'
#' @param sampleDirs Sample directories.
#'
#' @return [list].
.cellularBarcodesList <- function(sampleDirs) {
    files <- sampleDirs %>%
        file.path(paste(basename(.), "barcodes.tsv", sep = "-")) %>%
        setNames(names(sampleDirs))
    if (!all(file.exists(files))) {
        stop("Cellular barcode file missing", call. = FALSE)
    }
    message("Reading cellular barcode distributions")
    pblapply(seq_along(files), function(a) {
        .readCellularBarcodeFile(files[a])
    }) %>%
        setNames(names(sampleDirs))
}
