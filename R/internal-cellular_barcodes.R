#' Cellular barcode operations
#'
#' @rdname cellular_barcodes
#' @keywords internal
#'
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @param file Cellular barcode TSV file.
#' @param sample_dirs Sample directories.
#' @param list Cellular barcodes list.



#' @rdname cellular_barcodes
#' @return [tibble].
#' @details Modified version of Klein Lab MATLAB code.
.proportional_cb <- function(object) {
    metadata <- metadata(object)[["sample_metadata"]] %>%
        .[, meta_priority_cols]
    cellular_barcodes <- bcbio(object, "cellular_barcodes")
    lapply(seq_along(cellular_barcodes), function(a) {
        cb <- cellular_barcodes[[a]] %>%
            mutate(log10_reads = log10(.data[["reads"]]))
        cb_hist <- hist(cb[["log10_reads"]], n = 100L, plot = FALSE)
        # `fLog` in Klein Lab code
        counts <- cb_hist[["counts"]]
        # `xLog` in Klein Lab code
        mids <-  cb_hist[["mids"]]
        tibble(
            sample_id = names(cellular_barcodes)[[a]],
            log10_reads_per_cell = mids,
            proportion_of_cells = counts * (10L ^ mids) /
                sum(counts * (10L ^ mids)))
    }) %>%
        set_names(names(cellular_barcodes)) %>%
        bind_rows %>%
        left_join(metadata, by = "sample_id")
}



#' @rdname cellular_barcodes
#' @return [tibble].
.read_cellular_barcode_file <- function(file) {
    readFileByExtension(
        file,
        col_names = c("cellular_barcode", "reads"),
        col_types = "ci") %>%
        mutate(cellular_barcode = snake(.data[["cellular_barcode"]]))
}



#' @rdname cellular_barcodes
#' @return [list].
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
    }) %>% set_names(names(sample_dirs))
}



#' @rdname cellular_barcodes
#' @return [tibble].
.bind_cellular_barcodes <- function(list) {
    lapply(seq_along(list), function(a) {
        sample_id <- names(list)[[a]] %>% snake
        list[[a]] %>%
            mutate(sample_id = !!sample_id)
    }) %>%
        bind_rows %>%
        mutate(rowname = str_c(.data[["sample_id"]],
                               .data[["cellular_barcode"]],
                               sep = "_"))
}
