#' `bcbioSCDataSet` Object Accessors
#'
#' @rdname bcbio
#' @name bcbio
#' @author Michael Steinbaugh
#'
#' @param object Primary object.
#' @param type Type of count data to retrieve.
#'
#' @return [bcbioSCDataSet].
NULL



#' @rdname bcbio
#' @export
setMethod("bcbio", "bcbioSCDataSet", function(object, type) {
    if (type %in% names(slot(object, "callers"))) {
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
    function(object, type, value) {
        slot(object, "callers")[[type]] <- value
        validObject(object)
        object
    })
