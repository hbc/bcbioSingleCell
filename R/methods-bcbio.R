#' Additional bcbio Run Data Accessor
#'
#' @rdname bcbio
#' @name bcbio
#' @author Michael Steinbaugh
#'
#' @importFrom bcbioBase bcbio bcbio<-
#'
#' @inheritParams general
#'
#' @param type Type of count data to retrieve.
#'
#' @return [bcbioSingleCell].
#'
#' @examples
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
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
    signature("seurat"),
    function(object, type) {
        stopifnot(.hasSlot(object, "misc"))
        bcbio <- slot(object, "misc")[["bcbio"]]
        if (missing(type)) {
            return(bcbio)
        }
        if (type %in% names(bcbio)) {
            return(bcbio[[type]])
        } else {
            return(NULL)
        }
    }
)



# Assignment methods ===========================================================
#' @rdname bcbio
#' @export
setMethod(
    "bcbio<-",
    signature(
        object = "seurat",
        value = "ANY"
    ),
    function(object, type, value) {
        stopifnot(.hasSlot(object, "misc"))
        if (is.null(slot(object, "misc")[["bcbio"]])) {
            abort("seurat object was not generated with bcbioSingleCell")
        }
        slot(object, "misc")[["bcbio"]][[type]] <- value
        validObject(object)
        object
    })
