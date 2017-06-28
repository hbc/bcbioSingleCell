#' Cellular barcode operations
#'
#' @rdname cellular_barcodes
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param file Cellular barcode TSV file.
#' @param sample_dirs Sample directories.
.read_cellular_barcode_file <- function(file) {
    .read_file(file, col_names = c("cellular_barcode", "reads"), col_types = "ci")
}



#' @rdname cellular_barcodes
.cellular_barcodes <- function(sample_dirs) {
    files <- sample_dirs %>%
        file.path(paste(basename(.), "barcodes.tsv", sep = "-")) %>%
        set_names(names(sample_dirs))
    if (!all(file.exists(files))) {
        stop("Cellular barcode file missing")
    }
    message("Reading cellular barcode distributions")
    pblapply(seq_along(files), function(a) {
        .read_cellular_barcode_file(files[a])
    }
    ) %>% set_names(names(sample_dirs))
}



#' @rdname cellular_barcodes
#' @param list Cellular barcodes list
.bind_cellular_barcodes <- function(list) {
    lapply(seq_along(list), function(a) {
        sample_name <- names(list)[[a]]
        list[[a]] %>%
            mutate(rowname =
                       paste(sample_name,
                             .data[["cellular_barcode"]],
                             sep = ":"),
                   cellular_barcode = NULL)
    }
    ) %>%
        bind_rows
}
