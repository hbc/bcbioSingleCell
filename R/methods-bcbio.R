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
setMethod("bcbio", "bcbioSingleCellANY", function(object, type) {
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
    signature(object = "bcbioSingleCellANY", value = "ANY"),
    function(object, type, value) {
        slot(object, "bcbio")[[type]] <- value
        validObject(object)
        object
    })
