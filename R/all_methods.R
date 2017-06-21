#' Accessors for the sparse count matrix of a [bcbioSCDataSet] object
#'
#' @rdname bcbio
#' @docType methods
#' @name bcbio
#' @aliases bcbio bcbio,bcbioSCDataSet-method
#'   bcbio<-,bcbioSCDataSet,matrix-method
#'
#' @param object [bcbioSCDataSet] object.
#' @param value An integer matrix or other object.
#' @param type type of count data to retrieve
#' @param ... Matrix count data.
#'
#' @return Matrix/Object containing count data.
NULL



#' @rdname bcbio
#' @export
bcbio.bcbioSCDataSet <- function(object, type = "counts") {
    if (type == "counts")
        return(assays(object)[["counts"]])
    if (type %in% names(slot(object, "callers")))
        return(slot(object, "callers")[[type]])
    message(type, " not found")
}

#' @rdname bcbio
#' @export
setMethod("bcbio", signature(object = "bcbioSCDataSet"), bcbio.bcbioSCDataSet)

#' @rdname bcbio
#' @exportMethod "bcbio<-"
setReplaceMethod(
    "bcbio", signature(object = "bcbioSCDataSet", value = "ANY"),
    function(object, type = "counts", value) {
        if (type == "counts") {
            assays(object)[["counts"]] <- value
            validObject(object)
        } else {
            slot(object, "callers")[[type]] <- value
        }
        object
    })
