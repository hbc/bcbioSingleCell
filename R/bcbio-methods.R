#' Additional bcbio Run Data Accessor
#'
#' @name bcbio
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @inheritParams general
#' @param slot Slot name of data inside accessor.
#'
#' @return Various data types.
#'
#' @examples
#' # seurat ====
#' names(bcbio(seurat_small))
#'
#' # Assignment method support
#' bcbio(seurat_small, "metadata")[["stash"]] <- "XXX"
#' bcbio(seurat_small, "metadata")[["stash"]]
NULL



#' @importFrom bcbioBase bcbio
#' @export
bcbioBase::bcbio

#' @importFrom bcbioBase bcbio<-
#' @export
bcbioBase::`bcbio<-`



# Methods ======================================================================
#' @rdname bcbio
#' @export
setMethod(
    "bcbio",
    signature("seurat"),
    function(object, slot) {
        stopifnot(.hasSlot(object, "misc"))
        bcbio <- object@misc[["bcbio"]]
        if (missing(slot)) {
            return(bcbio)
        }
        if (!slot %in% names(bcbio)) {
            return(NULL)
        }
        bcbio[[slot]]
    }
)



#' @rdname bcbio
#' @export
setMethod(
    "bcbio<-",
    signature(
        object = "seurat",
        value = "ANY"
    ),
    function(object, slot, value) {
        stopifnot(.hasSlot(object, "misc"))
        if (is.null(object@misc[["bcbio"]])) {
            stop("seurat object was not generated with bcbioSingleCell")
        }
        object@misc[["bcbio"]][[slot]] <- value
        validObject(object)
        object
    }
)
