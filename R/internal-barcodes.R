#' Read sample cellular barcode summary file
#'
#' @rdname barcodes
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param file_name File name.
#' @param sample_dirs Sample directories.
.read_barcode_file <- function(file_name) {
    .read_file(file_name) %>%
        set_names(c("cellular_barcode", "reads"))
}

#' @rdname barcodes
.barcodes <- function(sample_dirs) {
    files <- file.path(sample_dirs,
                       paste(basename(.), "barcodes.tsv", sep = "-")) %>%
        set_names(names(sample_dirs))
    if (!all(file.exists(files))) {
        stop("Barcode TSV file missing")
    }
    message("Reading barcode distributions")
    pbmclapply(seq_along(files), function(a) {
        .read_barcode_file(files[a])
    }) %>% set_names(names(files))
}
