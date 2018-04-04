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
        data <- as.data.frame(rowData(object))
        assert_is_non_empty(data)
        cols <- c("geneID", "geneName")
        assert_is_subset(cols, colnames(data))
        data <- data[, cols]
        rownames(data) <- data[["geneID"]]
        data
    }
)



#' @rdname seurat
#' @importFrom bcbioBase interestingGroups
#' @export
#' @examples
#' # interestingGroups ====
#' interestingGroups(seurat_small)
#' interestingGroups(seurat_small) <- "sampleID"
#' interestingGroups(seurat_small)
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

#' @rdname seurat
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
#' @export
#' @examples
#' # rowData ====
#' rowData(pbmc_small) %>% glimpse()
#' rowData(seurat_small) %>% glimpse()
setMethod(
    "rowData",
    signature("seurat"),
    function(x) {
        rowRanges <- bcbio(x, "rowRanges")
        if (is(rowRanges, "GRanges")) {
            as(rowRanges, "DataFrame")
        } else {
            NULL
        }
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
