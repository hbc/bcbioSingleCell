#' Sample Barcode Metrics
#'
#' @rdname metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @return [data.frame] with cellular barcodes as rows.



#' @rdname metrics
#' @usage NULL
.metrics <- function(object) {
    meta <- sample_metadata(object) %>%
        tidy_select(unique(c(meta_priority_cols, interesting_groups(object))))
    col_data <- colData(object) %>%
        as.data.frame %>%
        rownames_to_column("cellular_barcode")
    # Sample identifier match fix
    if (metadata(object)[["pipeline"]] == "bcbio") {
        cb <- pull(col_data, "cellular_barcode")
        pattern <- "^(.*)_([acgt]{8})_([acgt]{8})_([acgt]{8})$"
        sample_ids <- str_match(cb, pattern) %>%
            as_tibble %>%
            set_colnames(c("cellular_barcode",
                           "file_name",
                           "revcomp",
                           "cb1",
                           "cb2")) %>%
            mutate(sample_id = paste(
                .data[["file_name"]],
                .data[["revcomp"]],
                sep = "_")) %>%
            pull("sample_id")
        col_data[["sample_id"]] <- sample_ids
    } else {
        col_data <- col_data %>%
            mutate(sample_id = .data[["cellular_barcode"]])
    }

    left_join(col_data, meta, by = "sample_id") %>%
        as.data.frame %>%
        column_to_rownames("cellular_barcode")
}



#' @rdname metrics
#' @export
setMethod("metrics", "bcbioSCDataSet", .metrics)



#' @rdname metrics
#' @export
setMethod("metrics", "bcbioSCSubset", .metrics)



#' @rdname metrics
#' @export
setMethod("metrics", "seurat", function(object) {
    object@data.info %>% snake
})
