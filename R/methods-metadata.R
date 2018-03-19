#' Metadata
#'
#' @name metadata
#'
#' @examples
#' # seurat ====
#' # Only works for objects created with bcbioSingleCell
#' names(metadata(seurat_small))
#'
#' # Assignment method support
#' metadata(seurat_small)[["stash"]] <- "XXX"
#' metadata(seurat_small)[["stash"]]
NULL



# Methods ======================================================================
#' @rdname metadata
#' @export
setMethod(
    "metadata",
    signature("seurat"),
    function(x) {
        bcbio(x, "metadata")
    }
)



# Assignment methods ===========================================================
#' @rdname metadata
#' @seealso `getMethod("metadata<-", "Annotated")`
#' @export
setMethod(
    "metadata<-",
    signature(
        x = "seurat",
        value = "ANY"
    ),
    function(x, value) {
        if (!is.list(value)) {
            abort("replacement 'metadata' value must be a list")
        }
        if (!length(value)) {
            names(value) <- NULL
        }
        bcbio(x, "metadata") <- value
        x
    }
)
