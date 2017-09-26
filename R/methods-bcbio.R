#' `bcbioSingleCell` Object Accessors
#'
#' @rdname bcbio
#' @name bcbio
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @param type Type of count data to retrieve.
#'
#' @return [bcbioSingleCell].
NULL



#' @rdname bcbio
#' @export
setMethod("bcbio", "bcbioSingleCell", function(object, type) {
    if (type %in% names(slot(object, "bcbio"))) {
        slot(object, "bcbio")[[type]]
    } else {
        stop(paste(type, "not found"))
    }
})



#' @rdname bcbio
#' @export
setMethod(
    "bcbio<-",
    signature(object = "bcbioSingleCell", value = "ANY"),
    function(object, type, value) {
        slot(object, "bcbio")[[type]] <- value
        validObject(object)
        object
    })



# Deprecated `bcbioSCDataSet` object support ====
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
