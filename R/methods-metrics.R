#' Sample barcode metrics
#'
#' @rdname metrics
#'
#' @param object Primary object.
#'
#' @return [data.frame].



#' @rdname metrics
#' @usage NULL
.metrics <- function(object) {
    meta <- sample_metadata(object) %>%
        tidy_select(unique(c(meta_priority_cols, interesting_groups(object))))
    colData(object) %>%
        as("data.frame") %>%
        rownames_to_column %>%
        mutate(
            sample_id = str_replace(.data[["rowname"]], "_[ACGT_]+$", "")) %>%
        left_join(meta, by = "sample_id") %>%
        column_to_rownames
}



#' @rdname metrics
#' @export
setMethod("metrics", "bcbioSCDataSet", .metrics)

#' @rdname metrics
#' @export
setMethod("metrics", "SCSubset", .metrics)

#' @rdname metrics
#' @export
setMethod("metrics", "seurat", function(object) {
    object@data.info %>% snake(rownames = FALSE)
})
