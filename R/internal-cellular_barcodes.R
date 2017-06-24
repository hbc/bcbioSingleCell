#' Read sample cellular barcode summary file
#'
#' @rdname barcodes
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param file_name File name.
#' @param sample_dirs Sample directories.
.read_cellular_barcode_file <- function(file_name) {
    .read_file(file_name) %>%
        set_names(c("cellular_barcode", "reads"))
}



#' Cellular barcode distributions
#'
#' @rdname cellular_barcodes
#' @keywords internal
#'
#' @author Michael Steinbaugh
.cellular_barcodes <- function(sample_dirs) {
    files <- sample_dirs %>%
        file.path(paste(basename(.), "barcodes.tsv", sep = "-")) %>%
        set_names(names(sample_dirs))
    if (!all(file.exists(files))) {
        stop("Cellular barcode file missing")
    }
    message("Reading cellular barcode distributions")
    pbmclapply(seq_along(files), function(a) {
        .read_cellular_barcode_file(files[a])
    }) %>% set_names(names(files))
}
