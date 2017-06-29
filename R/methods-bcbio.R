#' bcbioSCDataSet object accessors
#'
#' @rdname bcbio
#' @docType methods
#'
#' @param object Primary object.
#' @param type Type of count data to retrieve.
#' @param value An integer matrix or other object.
#'
#' @return Matrix/Object containing count data.
#' @export
setMethod("bcbio", "bcbioSCDataSet", function(object, type = "counts") {
    if (type == "counts") {
        assays(object)[["counts"]]
    } else if (type %in% names(slot(object, "callers"))) {
        slot(object, "callers")[[type]]
    } else {
        stop(paste(type, "not found"))
    }
})

#' @rdname bcbio
#' @export
setMethod(
    "bcbio<-",
    signature(object = "bcbioSCDataSet", value = "ANY"),
    function(object, type = "counts", value) {
        if (type == "counts") {
            assays(object)[["counts"]] <- value
            validObject(object)
        } else {
            slot(object, "callers")[[type]] <- value
        }
        object
    })
