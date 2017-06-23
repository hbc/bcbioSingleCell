#' Metadata table for knit report.
#'
#' @rdname metadata_table
#' @docType methods
#'
#' @author Michael Steinbaugh
#'
#' @param object Object.
#' @param ... Additional parameters.
#'
#' @return [kable].
#' @export
setMethod("metadata_table", "bcbioSCDataSet", function(object, ...) {
    colData(object) %>%
        as_tibble %>%
        remove_rownames %>%
        snake %>%
        mutate(file_name = NULL) %>%
        kable(caption = "Sample metadata")
})
