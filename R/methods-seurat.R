#' Extend S4 Methods for `seurat` Objects
#'
#' @name seurat
#' @author Michael Steinbaugh
#'
#' @examples
#' # dimensions ====
#' rownames(pbmc_small) %>% head()
#' colnames(pbmc_small) %>% head()
#'
#' # metadata ====
#' # Only works for objects created with bcbioSingleCell
#' names(metadata(seurat_small))
#'
#' # Assignment method support
#' metadata(seurat_small)[["stash"]] <- "XXX"
#' metadata(seurat_small)[["stash"]]
NULL



# Methods ======================================================================
#' @rdname seurat
#' @importFrom BiocGenerics colnames
#' @export
setMethod(
    "colnames",
    signature("seurat"),
    function(x) {
        colnames(counts(x))
    }
)



#' @rdname seurat
#' @export
setMethod(
    "metadata",
    signature("seurat"),
    function(x) {
        bcbio(x, "metadata")
    }
)



#' @rdname seurat
#' @importFrom BiocGenerics rownames
#' @export
setMethod(
    "rownames",
    signature("seurat"),
    function(x) {
        rownames(counts(x))
    }
)



# Assignment methods ===========================================================
#' @rdname seurat
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
