#' Sample barcode metrics
#'
#' @rdname metrics
#'
#' @param object [bcbioSCDataSet].
#'
#' @return [data.frame].
#' @export
setMethod("metrics", "bcbioSCDataSet", function(object) {
    interesting_groups <- interesting_groups(object)
    meta <- sample_metadata(object) %>%
        .[, unique(c("sample_name", interesting_groups))] %>%
        rownames_to_column("sample_id")
    colData(object) %>%
        as.data.frame %>%
        rownames_to_column %>%
        mutate(separate = .data[["rowname"]]) %>%
        separate_(col = "separate",
                  into = c("sample_id", "cellular_barcode"),
                  sep = ":") %>%
        left_join(meta, by = "sample_id") %>%
        mutate(cellular_barcode = NULL,
               sample_id = NULL) %>%
        column_to_rownames
})
