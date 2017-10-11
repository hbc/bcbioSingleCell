#' `bcbioSingleCell` Additional bcbio Run Data Accessor
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



# Methods ====
#' @rdname bcbio
#' @export
setMethod(
    "bcbio",
    signature("bcbioSingleCell"),
    function(object, type) {
        if (missing(type)) {
            return(slot(object, "bcbio"))
        }
        if (type %in% names(slot(object, "bcbio"))) {
            slot(object, "bcbio")[[type]]
        } else {
            stop(paste(type, "not found"), call. = FALSE)
        }
    })



#' @rdname bcbio
#' @export
setMethod(
    "bcbio<-",
    signature(object = "bcbioSingleCell",
              value = "ANY"),
    function(object, type, value) {
        slot(object, "bcbio")[[type]] <- value
        validObject(object)
        object
    })



# Legacy class support ====
# Package versions prior to 0.0.19 used `callers` to define the extra bcbio
# slot. The structure of the object is otherwise the same.
#' @rdname bcbio
#' @export
setMethod(
    "bcbio",
    signature("bcbioSCDataSet"),
    function(object, type) {
        if (missing(type)) {
            return(slot(object, "callers"))
        }
        if (type %in% names(slot(object, "callers"))) {
            slot(object, "callers")[[type]]
        } else {
            stop(paste(type, "not found"), call. = FALSE)
        }
    })



#' @rdname bcbio
#' @export
setMethod(
    "bcbio<-",
    signature(object = "bcbioSCDataSet",
              value = "ANY"),
    function(object, type, value) {
        slot(object, "callers")[[type]] <- value
        validObject(object)
        object
    })
