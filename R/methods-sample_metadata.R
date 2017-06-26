#' Sample metadata
#'
#' @rdname sample_metadata
#'
#' @param object Object.
#'
#' @export
setMethod("sample_metadata", "bcbioSCDataSet", function(object) {
    metadata(object)[["sample_metadata"]]
})
