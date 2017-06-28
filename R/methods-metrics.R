#' Sample barcode metrics
#'
#' @rdname metrics
#'
#' @param object [bcbioSCDataSet].
#' @param unique_names Unique sample names.
#'
#' @return [data.frame].



#' @rdname metrics
#' @usage NULL
.metrics <- function(object, unique_names = FALSE) {
    interesting_groups <- interesting_groups(object)
    meta <- sample_metadata(object, unique_names = unique_names) %>%
        tidy_select(unique(c("file_name",
                             "sample_id",
                             "sample_name",
                             interesting_groups)))
    metrics <- colData(object) %>%
        as.data.frame %>%
        rownames_to_column %>%
        mutate(separate = .data[["rowname"]]) %>%
        separate_(col = "separate",
                  into = c("sample_id", "cellular_barcode"),
                  sep = ":") %>%
        mutate(cellular_barcode = NULL)
    left_join(meta, metrics, by = "sample_id") %>%
        column_to_rownames
}



#' @rdname metrics
#' @export
setMethod("metrics", "bcbioSCDataSet", .metrics)

#' @rdname metrics
#' @export
setMethod("metrics", "SummarizedExperiment", .metrics)
