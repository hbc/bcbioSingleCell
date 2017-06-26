#' Sample metadata
#'
#' @rdname sample_metadata
#'
#' @author Michael Steinbaugh
#'
#' @param object Object.
#'
#' @export
setMethod("sample_metadata", "bcbioSCDataSet", function(object) {
    metadata(object)[["sample_metadata"]] %>% as.data.frame
})

#' @rdname sample_metadata
#' @export
setMethod("sample_metadata", "bcbioSCSubset", function(object) {
    object@callers[["sample_metadata"]] %>% as.data.frame
})
