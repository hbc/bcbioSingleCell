#' Additional bcbio Run Data Accessor
#'
#' @rdname bcbio
#' @name bcbio
#' @author Michael Steinbaugh
#'
#' @importFrom basejump bcbio bcbio<-
#'
#' @inheritParams AllGenerics
#'
#' @param type Type of count data to retrieve.
#'
#' @return [bcbioSingleCell].
#'
#' @examples
#' load(system.file(
#'     file.path("inst", "extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("inst", "extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' bcbio(bcb) %>% names()
#'
#' # seurat
#' bcbio(seurat) %>% names()
NULL



# Methods ====
#' @rdname bcbio
#' @export
setMethod(
    "bcbio",
    signature("bcbioSingleCell"),
    function(object, type) {
        bcbio <- slot(object, "bcbio")
        if (missing(type)) {
            return(bcbio)
        }
        if (type %in% names(bcbio)) {
            return(bcbio[[type]])
        } else {
            return(NULL)
        }
    })



#' @rdname bcbio
#' @export
setMethod(
    "bcbio",
    signature("seurat"),
    function(object, type) {
        if (!.hasSlot(object, "misc")) {
            return(NULL)
        }
        bcbio <- slot(object, "misc")[["bcbio"]]
        if (missing(type)) {
            return(bcbio)
        }
        if (type %in% names(bcbio)) {
            return(bcbio[[type]])
        } else {
            return(NULL)
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
            return(slot(object, "callers")[[type]])
        } else {
            return(NULL)
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
