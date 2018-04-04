#' Extend S4 Methods for `seurat` Objects
#'
#' Provide `SummarizedExperiment` method support.
#'
#' @name seurat
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return Match `SummarizedExperiment` method return.
NULL



# Methods ======================================================================
#' @rdname seurat
#' @importFrom SummarizedExperiment assay
#' @export
#' @examples
#' # assay ====
#' assay(pbmc_small) %>% glimpse()
setMethod(
    "assay",
    signature("seurat"),
    function(x) {
        slot(x, "raw.data")
    }
)



#' @rdname seurat
#' @importFrom BiocGenerics colnames
#' @export
#' @examples
#' # colnames ====
#' colnames(pbmc_small) %>% head()
setMethod(
    "colnames",
    signature("seurat"),
    function(x) {
        colnames(counts(x))
    }
)



#' @rdname seurat
#' @importFrom bcbioBase gene2symbol
#' @export
#' @examples
#' # gene2symbol ====
#' gene2symbol(seurat_small) %>% glimpse()
setMethod(
    "gene2symbol",
    signature("seurat"),
    function(object) {
        rowData <- rowData(object)
        assert_is_non_empty(rowData)
        cols <- c("geneID", "geneName")
        assert_is_subset(cols, colnames(rowData))
        rowData[, cols]
    }
)



#' @rdname seurat
#' @export
#' @examples
#' # metadata ====
#' names(metadata(seurat_small))
#' metadata(seurat_small)[["stash"]] <- "XXX"
#' metadata(seurat_small)[["stash"]]
setMethod(
    "metadata",
    signature("seurat"),
    function(x) {
        bcbio(x, "metadata")
    }
)



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



#' @rdname seurat
#' @importFrom BiocGenerics rownames
#' @export
#' @examples
#' # rownames ====
#' rownames(pbmc_small) %>% head()
setMethod(
    "rownames",
    signature("seurat"),
    function(x) {
        rownames(counts(x))
    }
)
