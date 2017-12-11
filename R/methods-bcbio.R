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
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' bcbio(bcb) %>% names()
#'
#' # Assignment method support
#' bcbio(bcb, "stash") <- "testing"
#' bcbio(bcb, "stash")
#'
#' # seurat
#' bcbio(seurat) %>% names()
#'
#' # Assignment method support
#' bcbio(seurat, "stash") <- "testing"
#' bcbio(seurat, "stash")
NULL



# Methods ======================================================================
#' @rdname bcbio
#' @export
setMethod(
    "bcbio",
    signature("bcbioSingleCell"),
    function(object, type) {
        bcbio <- slot(object, "bcbio")
        if (missing(type)) return(bcbio)
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
        if (!.hasSlot(object, "misc")) return(NULL)
        bcbio <- slot(object, "misc")[["bcbio"]]
        if (missing(type)) return(bcbio)
        if (type %in% names(bcbio)) {
            return(bcbio[[type]])
        } else {
            return(NULL)
        }
    })



# Assignment methods ===========================================================
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



#' @rdname bcbio
#' @export
setMethod(
    "bcbio<-",
    signature(object = "seurat",
              value = "ANY"),
    function(object, type, value) {
        if (!.hasSlot(object, "misc")) return(NULL)
        if (is.null(slot(object, "misc")[["bcbio"]])) {
            stop("seurat object was not generated with bcbioSingleCell",
                 call. = FALSE)
        }
        slot(object, "misc")[["bcbio"]][[type]] <- value
        validObject(object)
        object
    })



# Legacy methods ===============================================================
# Package versions prior to 0.0.19 used `callers` to define the extra bcbio
# slot. The structure of the object is otherwise the same.
#' @rdname bcbio
#' @export
setMethod(
    "bcbio",
    signature("bcbioSCDataSet"),
    function(object, type) {
        if (missing(type)) return(slot(object, "callers"))
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
