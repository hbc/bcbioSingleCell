#' S4 generics
#'
#' @rdname generics
#'
#' @param object Object.
#' @param value Value to assign.
#' @param ... Additional parameters.



#' @rdname generics
#' @export
setGeneric("aggregate_replicates", function(object, ...) {
    standardGeneric("aggregate_replicates")
})



#' @rdname generics
#' @export
setGeneric("bcbio", function(object, ...) standardGeneric("bcbio"))



#' @rdname generics
#' @export
setGeneric("bcbio<-", function(object, ..., value) standardGeneric("bcbio<-"))



#' @rdname generics
#' @export
setGeneric("filter_barcodes", function(object, ...) {
    standardGeneric("filter_barcodes")
})



#' @rdname generics
#' @export
setGeneric("metadata_table", function(object, ...) {
    standardGeneric("metadata_table")
})



#' @rdname generics
#' @export
setGeneric("interesting_groups", function(object) {
    standardGeneric("interesting_groups")
})



#' @rdname generics
#' @export
setGeneric("metrics", function(object, ...) standardGeneric("metrics"))



#' @rdname generics
#' @export
setGeneric("sample_metadata", function(object, ...) {
    standardGeneric("sample_metadata")
})



#' @rdname generics
#' @export
setGeneric("select_samples", function(object, ...) {
    standardGeneric("select_samples")
})



#' @rdname generics
#' @export
setGeneric("top_barcodes", function(object, ...) {
    standardGeneric("top_barcodes")
})
