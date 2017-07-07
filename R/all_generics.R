#' @rdname aggregate_replicates
#' @usage NULL
#' @export
setGeneric("aggregate_replicates", function(object, ...) {
    standardGeneric("aggregate_replicates")
})



#' @rdname bcbio
#' @usage NULL
#' @export
setGeneric("bcbio", function(object, ...) standardGeneric("bcbio"))



#' @rdname bcbio
#' @usage NULL
#' @export
setGeneric("bcbio<-", function(object, ..., value) standardGeneric("bcbio<-"))



#' @rdname metadata_table
#' @usage NULL
#' @export
setGeneric("metadata_table", function(object, ...) {
    standardGeneric("metadata_table")
})



#' @rdname interesting_groups
#' @usage NULL
#' @export
setGeneric("interesting_groups", function(object) {
    standardGeneric("interesting_groups")
})



#' @rdname metrics
#' @usage NULL
#' @export
setGeneric("metrics", function(object, ...) standardGeneric("metrics"))



#' @rdname sample_metadata
#' @usage NULL
#' @export
setGeneric("sample_metadata", function(object, ...) {
    standardGeneric("sample_metadata")
})



#' @rdname select_samples
#' @usage NULL
#' @export
setGeneric("select_samples", function(object, ...) {
    standardGeneric("select_samples")
})



#' @rdname top_barcodes
#' @usage NULL
#' @export
setGeneric("top_barcodes", function(object, ...) {
    standardGeneric("top_barcodes")
})
