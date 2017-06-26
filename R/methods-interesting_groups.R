#' Interesting groups
#'
#' @rdname interesting_groups
#'
#' @param object Object.
#'
#' @export
setMethod("interesting_groups", "bcbioSCDataSet", function(object) {
    metadata(object)[["interesting_groups"]]
})



#' @rdname interesting_groups
#' @export
setMethod("interesting_groups", "bcbioSCSubset", function(object) {
    object@callers[["interesting_groups"]]
})
