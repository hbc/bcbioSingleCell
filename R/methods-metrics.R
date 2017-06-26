#' Sample barcode metrics
#'
#' @rdname metrics
#'
#' @param object [bcbioSCDataSet].
#'
#' @return [data.frame].
#' @export
setMethod("metrics", "bcbioSCDataSet", function(object) {
    interesting_groups <- metadata(object)[["interesting_groups"]]
    metadata <- sample_metadata(object) %>%
        .[, unique(c("sample_name", interesting_groups))]
    metrics <- colData(object) %>%
        as.data.frame %>%
        rownames_to_column %>%
        mutate(separate = .data[["rowname"]]) %>%
        separate_(col = "separate",
                  into = c("sample_name", "barcode"),
                  sep = "-") %>%
        mutate(barcode = NULL) %>%
        left_join(metadata, by = "sample_name") %>%
        column_to_rownames
})
