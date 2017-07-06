#' Sample barcode metrics
#'
#' @rdname metrics
#'
#' @param object Primary object.
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
    colData(object) %>%
        as.data.frame %>%
        rownames_to_column %>%
        mutate(sample_id = str_replace(.data[["rowname"]],
                                       "_[acgt_]+$",
                                       "")) %>%
    left_join(meta, by = "sample_id") %>%
        column_to_rownames
}



#' @rdname metrics
#' @export
setMethod("metrics", "bcbioSCDataSet", .metrics)

#' @rdname metrics
#' @export
setMethod("metrics", "SummarizedExperiment", .metrics)
