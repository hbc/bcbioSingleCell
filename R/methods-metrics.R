#' Sample Barcode Metrics
#'
#' @rdname metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @return [data.frame] with cellular barcodes as rows.



#' @rdname metrics
#' @usage NULL
.metrics <- function(object) {
    umi_type <- metadata(object)[["umi_type"]]
    meta <- sample_metadata(object) %>%
        tidy_select(unique(c(meta_priority_cols, interesting_groups(object))))
    col_data <- colData(object) %>%
        as.data.frame %>%
        rownames_to_column
    if (umi_type == "surecell") {
        match <- str_match(col_data[["rowname"]],
                           "^(.+)_([acgt]{6}_[acgt]{6}_[acgt]{6})$")
    } else {
        match <- str_match(col_data[["rowname"]],
                           "^(.+)_([acgt]{8}_[acgt]{8})$")
    }
    match <- match %>%
        as.data.frame %>%
        set_colnames(c("rowname", "sample_id", "cellular_barcode"))

    left_join(col_data, match, by = "rowname") %>%
        left_join(meta, by = "sample_id") %>%
        as.data.frame %>%
        column_to_rownames
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
