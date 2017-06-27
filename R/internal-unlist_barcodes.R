#' Unlist cellular barcodes
#'
#' Convert named list of cellular barcodes per sample to a data frame.
#'
#' @rdname unlist_barcodes
#' @keywords internal
#'
#' @author Michael Steinbaugh
#'
#' @param bcb [bcbioSCDataSet].
.unlist_barcodes <- function(bcb) {
    barcodes <- bcbio(bcb, "cellular_barcodes")
    meta <- sample_metadata(bcb) %>%
        .[, c("sample_id", "sample_name")]
    message("Converting nested barcodes to data frame")
    pbmclapply(seq_along(barcodes), function(a) {
        barcodes[a] %>%
            as.data.frame %>%
            rownames_to_column %>%
            set_names(c("cellular_barcode", "reads")) %>%
            arrange(!!!syms(c("reads", "cellular_barcode"))) %>%
            mutate(log10_reads = log10(.data[["reads"]]),
                   sample_id = names(barcodes[a]))
    }
    ) %>%
        bind_rows %>%
        left_join(meta, by = "sample_id")
}
