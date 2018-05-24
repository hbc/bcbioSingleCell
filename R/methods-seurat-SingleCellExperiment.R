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
        x <- as(x, "SingleCellExperiment")
        assay(x)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @export
setMethod(
    "colData",
    signature("seurat"),
    function(x) {
        x <- as(x, "SingleCellExperiment")
        colData(x)
    }
)



# Note that Seurat subset operations keep `raw.data` matrix unmodified by
# default and only subset the `data` matrix
#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics colnames
#' @export
setMethod(
    "colnames",
    signature("seurat"),
    function(x) {
        x <- as(x, "SingleCellExperiment")
        colnames(x)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics counts
#' @param normalized Normalized (`TRUE`) or raw (`FALSE`) counts.
#' @export
setMethod(
    "counts",
    signature("seurat"),
    function(object) {
        x <- as(x, "SingleCellExperiment")
        counts(x)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "metadata",
    signature("seurat"),
    function(x) {
        stash <- bcbio(x, "metadata")
        if (!is.null(stash)) {
            message("Using bcbio stashed metadata")
            return(stash)
        }
        x <- as(x, "SingleCellExperiment")
        metadata(x)
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
        rr <- rowRanges(x)
        x <- as(x, "SingleCellExperiment")
        rowRanges(x) <- rr
        rowData(x)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "rowRanges",
    signature("seurat"),
    function(x) {
        stash <- bcbio(x, "rowRanges")
        x <- as(x, "SingleCellExperiment")
        if (is(stash, "GRanges")) {
            assert_are_identical(names(stash), rownames(x))
            message("Using bcbio stashed rowRanges")
            return(stash)
        }
        x <- as(x, "SingleCellExperiment")
        rowRanges(x)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics rownames
#' @export
setMethod(
    "rownames",
    signature("seurat"),
    function(x) {
        x <- as(x, "SingleCellExperiment")
        rownames(x)
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



#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioBase sampleNames
#' @export
setMethod(
    "sampleNames",
    signature("seurat"),
    getMethod("sampleNames", "SummarizedExperiment")
)
