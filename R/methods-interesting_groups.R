#' Interesting Groups
#'
#' @rdname interesting_groups
#'
#' @return Character vector.



#' @rdname interesting_groups
#' @usage NULL
.interesting_groups <- function(object) {
    metadata(object)[["interesting_groups"]]
}



#' @rdname interesting_groups
#' @export
setMethod("interesting_groups", "bcbioSCDataSet", .interesting_groups)



#' @rdname interesting_groups
#' @export
setMethod("interesting_groups", "bcbioSCSubset", .interesting_groups)
