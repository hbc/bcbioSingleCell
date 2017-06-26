#' @rdname aggregate_replicates
#' @export
setGeneric("aggregate_replicates", function(object, ...) {
    standardGeneric("aggregate_replicates")
})



#' @rdname bcbio
#' @export
setGeneric("bcbio", function(object, ...) standardGeneric("bcbio"))



#' @rdname bcbio
#' @export
setGeneric("bcbio<-", function(object, ..., value) standardGeneric("bcbio<-"))



#' @rdname counts
counts <- BiocGenerics::counts



#' @rdname filter_barcodes
#' @export
setGeneric("filter_barcodes", function(object, ...) {
    standardGeneric("filter_barcodes")
})



#' @rdname metadata_table
#' @export
setGeneric("metadata_table", function(object, ...) {
    standardGeneric("metadata_table")
})



#' @rdname interesting_groups
#' @export
setGeneric("interesting_groups", function(object) {
    standardGeneric("interesting_groups")
})



#' @rdname metrics
#' @export
setGeneric("metrics", function(object, ...) standardGeneric("metrics"))



#' @rdname sample_metadata
#' @export
setGeneric("sample_metadata", function(object) {
    standardGeneric("sample_metadata")
})



#' @rdname select_samples
#' @export
setGeneric("select_samples", function(object, ...) {
    standardGeneric("select_samples")
})



#' @rdname top_barcodes
#' @export
setGeneric("top_barcodes", function(object, ...) {
    standardGeneric("top_barcodes")
})
