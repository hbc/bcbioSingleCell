#' Extend S4 Methods for `seurat` Class
#'
#' Provide limited `SingleCellExperiment`-like method support.
#'
#' @name seurat-SingleCellExperiment
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return Match `SummarizedExperiment` method return.
NULL



# SingleCellExperiment =========================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @export
setMethod(
    "assay",
    signature("seurat"),
    function(x) {
        slot(x, "raw.data")
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @export
setMethod(
    "colData",
    signature("seurat"),
    function(x) {
        data <- slot(x, "meta.data")
        as(data, "DataFrame")
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics colnames
#' @export
setMethod(
    "colnames",
    signature("seurat"),
    function(x) {
        colnames(counts(x))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "metadata",
    signature("seurat"),
    function(x) {
        bcbio(x, "metadata")
    }
)



#' @rdname seurat-SingleCellExperiment
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
            stop("replacement 'metadata' value must be a list")
        }
        if (!length(value)) {
            names(value) <- NULL
        }
        bcbio(x, "metadata") <- value
        x
    }
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "rowData",
    signature("seurat"),
    function(x) {
        rowRanges <- rowRanges(x)
        if (is.null(rowRanges)) {
            return(NULL)
        }
        as(rowRanges, "DataFrame")
    }
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "rowRanges",
    signature("seurat"),
    function(x) {
        bcbio(x, "rowRanges")
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics rownames
#' @export
setMethod(
    "rownames",
    signature("seurat"),
    function(x) {
        rownames(counts(x))
    }
)



# bcbioBase Methods ============================================================
#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioBase gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("seurat"),
    getMethod("gene2symbol", "SummarizedExperiment")
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioBase interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("seurat"),
    function(object) {
        validObject(object)
        value <- metadata(object)[["interestingGroups"]]
        if (is.null(value)) {
            value <- "sampleName"
        }
        assertFormalInterestingGroups(object, value)
        value
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioBase interestingGroups<-
#' @export
setMethod(
    "interestingGroups<-",
    signature(
        object = "seurat",
        value = "character"
    ),
    getMethod(
        "interestingGroups<-",
        signature(
            object = "SummarizedExperiment",
            value = "character"
        )
    )
)
