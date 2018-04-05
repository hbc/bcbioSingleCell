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
            abort("replacement 'metadata' value must be a list")
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
        if (is(rowRanges, "GRanges")) {
            as(rowRanges, "DataFrame")
        } else {
            NULL
        }
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



# bcbioSingleCell Methods ======================================================
#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "cell2sample",
    signature("seurat"),
    getMethod("cell2sample", "SingleCellExperiment")
)



#' @rdname diffExp
#' @export
setMethod(
    "diffExp",
    signature("seurat"),
    getMethod("diffExp", "SingleCellExperiment")
)



#' @rdname fetchGeneData
#' @export
setMethod(
    "fetchGeneData",
    signature("seurat"),
    getMethod("fetchGeneData", "SingleCellExperiment")
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioBase gene2symbol
#' @export
setMethod(
    "gene2symbol",
    signature("seurat"),
    function(object) {
        data <- as.data.frame(rowData(object))
        assert_is_non_empty(data)
        cols <- c("geneID", "geneName")
        assert_is_subset(cols, colnames(data))
        data <- data[, cols]
        rownames(data) <- data[["geneID"]]
        data
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom bcbioBase interestingGroups
#' @export
setMethod(
    "interestingGroups",
    signature("seurat"),
    function(object) {
        validObject(object)
        x <- metadata(object)[["interestingGroups"]]
        if (is.null(x)) {
            x <- "sampleName"
        }
        x
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
    function(object, value) {
        assertFormalInterestingGroups(
            x = sampleData(object),
            interestingGroups = value
        )
        if (is.null(metadata(object))) {
            abort("object was not created with bcbioSingleCell")
        }
        metadata(object)[["interestingGroups"]] <- value
        validObject(object)
        object
    }
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "plotCellCounts",
    signature("seurat"),
    getMethod("plotCellCounts", "SingleCellExperiment")
)



#' @rdname plotCumulativeUMIsPerCell
#' @export
setMethod(
    "plotCumulativeUMIsPerCell",
    signature("seurat"),
    getMethod("plotCumulativeUMIsPerCell", "SingleCellExperiment")
)



#' @rdname plotReadsPerCell
#' @export
setMethod(
    "plotReadsPerCell",
    signature("seurat"),
    getMethod("plotReadsPerCell", "SingleCellExperiment")
)

#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("seurat"),
    getMethod("plotZerosVsDepth", "SingleCellExperiment")
)
