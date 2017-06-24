#' Unlist cellular barcodes
#'
#' Convert named list of cellular barcodes per sample to a data frame.
#'
#' @rdname unlist_barcodes
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param run [bcbioSCDataSet].
.unlist_barcodes <- function(run) {
    barcodes <- run$barcodes
    message("Converting nested barcodes to data frame...")
    pbmclapply(seq_along(barcodes), function(a) {
        barcodes[a] %>%
            as.data.frame %>%
            rownames_to_column %>%
            set_names(c("cellular_barcode", "reads")) %>%
            arrange(!!!syms(c("reads", "cellular_barcode"))) %>%
            mutate(log10_reads = log10(.data$reads),
                   sample_barcode = names(barcodes[a]))
    }) %>% bind_rows %>%
        left_join(run$metadata[, c("sample_barcode", "sample_name")],
                  by = "sample_barcode")
}
