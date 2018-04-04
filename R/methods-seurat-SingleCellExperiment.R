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



# Methods ======================================================================
#' @rdname seurat-SingleCellExperiment
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



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @export
#' @examples
#' # colData ====
#' colData(seurat_small) %>% glimpse()
#' colData(pbmc_small) %>% glimpse()
setMethod(
    "colData",
    signature("seurat"),
    function(x) {
        sampleData <- sampleData(x, return = "DataFrame")
        colData <- slot(x, "meta.data")
        assert_is_data.frame(colData)
        colData <- as(colData, "DataFrame")
        assert_are_disjoint_sets(
            x = colnames(sampleData),
            y = colnames(colData)
        )
        colData <- camel(colData)
        cell2sample <- cell2sample(x)
        sampleID <- DataFrame("sampleID" = cell2sample)
        colData <- cbind(colData, sampleID)
        colData[["rowname"]] <- rownames(colData)
        data <- merge(
            x = colData,
            y = sampleData,
            by = "sampleID",
            all.x = TRUE
        )
        # Ensure the numeric metrics columns appear first
        data <- data[, unique(c(colnames(colData), colnames(data)))]
        rownames(data) <- data[["rowname"]]
        data[["rowname"]] <- NULL
        # Add `interestingGroups` column
        interestingGroups <- interestingGroups(x)
        data <- uniteInterestingGroups(data, interestingGroups)
        # Add `ident` column
        ident <- slot(x, "ident")
        assert_is_factor(ident)
        ident <- DataFrame("ident" = ident)
        data <- cbind(data, ident)
        data
    }
)



#' @rdname seurat-SingleCellExperiment
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



#' @rdname seurat-SingleCellExperiment
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



#' @rdname seurat-SingleCellExperiment
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
#' @examples
#' # rowData ====
#' rowData(seurat_small) %>% glimpse()
#' rowData(pbmc_small) %>% glimpse()
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
#' @examples
#' # rowRanges ====
#' rowRanges(seurat_small) %>% glimpse()
#' rowRanges(pbmc_small) %>% glimpse()
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
