#' @rdname bcbio
#' @export
setGeneric("bcbio", function(object, ...) standardGeneric("bcbio"))

#' @rdname bcbio
#' @export
setGeneric("bcbio<-", function(object, ..., value) standardGeneric("bcbio<-"))

#' @rdname metadata_table
#' @export
setGeneric("metadata_table", function(object, ...) {
    standardGeneric("metadata_table")
})
