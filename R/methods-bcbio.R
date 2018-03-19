# FIXME Update handling of gene2symbol and rowRanges here

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
#' # seurat ====
#' names(bcbio(seurat_small))
#'
#' # Assignment method support
#' bcbio(seurat_small, "stash") <- "testing"
#' bcbio(seurat_small, "stash")
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
    }
)
