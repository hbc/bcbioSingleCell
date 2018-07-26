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



# Assert check to see if we're modifying a freshly created seurat object
.isNewSeurat <- function(object) {
    assert_are_identical(x@raw.data, x@data)
    stopifnot(is.null(x@scale.data))
    stopifnot(!length(x@var.genes))
}



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



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment colData<-
#' @export
setMethod(
    "colData<-",
    signature(
        x = "seurat",
        value = "DataFrame"
    ),
    function(x, value) {
        slot(x, "meta.data") <- as.data.frame(value)
        validObject(x)
        x
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
#' @importFrom BiocGenerics colnames<-
#' @export
setMethod(
    "colnames<-",
    signature("seurat"),
    function(x, value) {
        .isNewSeurat(x)
        colnames(x@raw.data) <- value
        x@data <- x@raw.data
        x
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump convertGenesToSymbols
#' @export
setMethod(
    "convertGenesToSymbols",
    signature("seurat"),
    getMethod("convertGenesToSymbols", "SummarizedExperiment")
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics counts
#' @param normalized Normalized (`TRUE`) or raw (`FALSE`) counts.
#' @export
setMethod(
    "counts",
    signature("seurat"),
    function(object) {
        object <- as(object, "SingleCellExperiment")
        counts(object)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("seurat"),
    getMethod("gene2symbol", "SummarizedExperiment")
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("seurat"),
    getMethod("interestingGroups", "SummarizedExperiment")
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump interestingGroups<-
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
#' @export
setMethod(
    "metadata",
    signature("seurat"),
    function(x) {
        stash <- slot(x, "misc")[["bcbio"]][["metadata"]]
        if (!is.null(stash)) {
            stash
        } else {
            x %>%
                as.SingleCellExperiment() %>%
                metadata()
        }
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
        slot(x, "misc")[["bcbio"]][["metadata"]] <- value
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
        stash <- slot(x, "misc")[["bcbio"]][["rowRanges"]]
        x <- as(x, "SingleCellExperiment")
        if (is(stash, "GRanges")) {
            assert_are_identical(names(stash), rownames(x))
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



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics rownames<-
#' @export
setMethod(
    "rownames<-",
    signature("seurat"),
    function(x, value) {
        .isNewSeurat(x)
        rownames(x@raw.data) <- value
        x@data <- x@raw.data
        x
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump sampleNames
#' @export
setMethod(
    "sampleNames",
    signature("seurat"),
    getMethod("sampleNames", "SummarizedExperiment")
)
