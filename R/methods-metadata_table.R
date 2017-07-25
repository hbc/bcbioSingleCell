#' Metadata Table
#'
#' Table output for a knit report.
#'
#' @rdname metadata_table
#' @author Michael Steinbaugh
#'
#' @param ... Additional arguments, passed to [knitr::kable()].
#'
#' @return [kable].



#' @rdname metadata_table
#' @usage NULL
.metadata_table <- function(object, ...) {
    object %>%
        sample_metadata %>%
        mutate(file_name = NULL) %>%
        kable(caption = "Sample metadata", ...)
}



#' @rdname metadata_table
#' @export
setMethod("metadata_table", "bcbioSCDataSet", .metadata_table)



#' @rdname metadata_table
#' @export
setMethod("metadata_table", "bcbioSCSubset", .metadata_table)
