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
