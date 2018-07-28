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
.assertIsNewSeurat <- function(object) {
    assert_are_identical(object@raw.data, object@data)
    stopifnot(is.null(object@scale.data))
    stopifnot(!length(object@var.genes))
}



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assay
#' @export
setMethod(
    "assay",
    signature("seurat"),
    function(x) {
        assay(as.SingleCellExperiment(x))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assays
#' @export
setMethod(
    "assays",
    signature("seurat"),
    function(x) {
        assays(as.SingleCellExperiment(x))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @export
setMethod(
    "colData",
    signature("seurat"),
    function(x) {
        colData(as.SingleCellExperiment(x))
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



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics colnames
#' @export
setMethod(
    "colnames",
    signature("seurat"),
    function(x) {
        colnames(as.SingleCellExperiment(x))
    }
)



#' Can't inherit SummarizedExperiment because seurat doesn't support dim, so
#' therefore we can't assign rownames
#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump convertGenesToSymbols
#' @export
setMethod(
    "convertGenesToSymbols",
    signature("seurat"),
    function(object) {
        validObject(object)
        .assertIsNewSeurat(object)
        gene2symbol <- gene2symbol(object)
        if (is.null(gene2symbol)) {
            warning("Object doesn't contain gene-to-symbol mappings")
            return(object)
        }
        symbols <- gene2symbol %>%
            .[, "geneName", drop = TRUE] %>%
            as.character() %>%
            make.unique()
        rownames(object@raw.data) <- symbols
        object@data <- object@raw.data
        object
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
        counts(as.SingleCellExperiment(x))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics colnames
#' @export
setMethod(
    "colnames",
    signature("seurat"),
    function(x) {
        colnames(as.SingleCellExperiment(x))
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
#' @importFrom SingleCellExperiment reducedDims
#' @export
setMethod(
    "reducedDims",
    signature("seurat"),
    function(x) {
        reducedDims(as.SingleCellExperiment(x))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "rowData",
    signature("seurat"),
    function(x) {
        rowData(as.SingleCellExperiment(x))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "rowRanges",
    signature("seurat"),
    function(x) {
        # Attempt to use stashed rowRanges, if present
        gr <- slot(x, "misc")[["bcbio"]][["rowRanges"]]
        if (is(gr, "GRanges")) {
            assert_is_subset(c("geneID", "geneName"), colnames(mcols(gr)))
            names(gr) <- make.unique(as.character(mcols(gr)[["geneName"]]))
            assert_is_subset(rownames(x), names(gr))
            return(gr[rownames(x)])
        }
        rowRanges(as.SingleCellExperiment(x))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics rownames
#' @export
setMethod(
    "rownames",
    signature("seurat"),
    function(x) {
        rownames(as.SingleCellExperiment(x))
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
