#' Metadata table for knit report
#'
#' @rdname metadata_table
#' @docType methods
#'
#' @author Michael Steinbaugh
#'
#' @param object Primary object.
#' @param ... Passthrough parameters to [knitr::kable()].
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
setMethod("metadata_table", "SCSubset", .metadata_table)
