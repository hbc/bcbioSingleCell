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



# bcbio Methods ================================================================
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



#' @rdname inflectionPoint
#' @export
setMethod(
    "inflectionPoint",
    signature("seurat"),
    getMethod("inflectionPoint", "SingleCellExperiment")
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
        if (is.character(x)) {
            x
        } else {
            NULL
        }
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



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("seurat"),
    function(object, ...) {
        fun <- getMethod("metrics", "SingleCellExperiment")
        data <- fun(object, ...)
        # Add ident column
        data[["ident"]] <- slot(object, "ident")
        data
    }
)



#' @rdname metricsPerSample
#' @export
setMethod(
    "metricsPerSample",
    signature("seurat"),
    getMethod("metricsPerSample", "SingleCellExperiment")
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



#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("seurat"),
    getMethod("plotGenesPerCell", "SingleCellExperiment")
)



#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("seurat"),
    getMethod("plotMitoRatio", "SingleCellExperiment")
)



#' @rdname plotMitoVsCoding
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("seurat"),
    getMethod("plotMitoVsCoding", "SingleCellExperiment")
)



#' @rdname plotQC
#' @export
setMethod(
    "plotQC",
    signature("seurat"),
    getMethod("plotQC", "SingleCellExperiment")
)



#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    signature("seurat"),
    getMethod("plotNovelty", "SingleCellExperiment")
)



#' @rdname plotQuantileHeatmap
#' @export
setMethod(
    "plotQuantileHeatmap",
    signature("seurat"),
    getMethod("plotQuantileHeatmap", "SummarizedExperiment")
)



#' @rdname plotReadsPerCell
#' @export
setMethod(
    "plotReadsPerCell",
    signature("seurat"),
    getMethod("plotReadsPerCell", "SingleCellExperiment")
)



#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("seurat"),
    getMethod("plotUMIsPerCell", "SingleCellExperiment")
)



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("seurat"),
    getMethod("plotUMIsVsGenes", "SingleCellExperiment")
)



#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("seurat"),
    getMethod("plotZerosVsDepth", "SingleCellExperiment")
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("seurat"),
    function(
        object,
        interestingGroups,
        return = c("DataFrame", "data.frame")
    ) {
        validObject(object)
        stopifnot(.hasSlot(object, "version"))
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        if (!is.character(interestingGroups)) {
            interestingGroups <- "sampleName"
        }
        return <- match.arg(return)
        data <- metadata(object)[["sampleData"]]
        if (is.null(data)) {
            data <- slot(object, "meta.data")
            assert_is_data.frame(data)
            # Create priority columns from `orig.ident`, if necessary
            if (!all(bcbioBase::metadataPriorityCols %in% colnames(data))) {
                missing <- setdiff(
                    x = bcbioBase::metadataPriorityCols,
                    y = colnames(data)
                )
                for (i in seq_along(missing)) {
                    data[[missing[[i]]]] <- data[["orig.ident"]]
                }
            }
            blacklist <- paste(
                c(
                    "cellularBarcode",
                    "orig.ident",
                    "Phase",
                    "^res\\.[0-9]"
                ),
                collapse = "|"
            )
            data <- data %>%
                remove_rownames() %>%
                .[, !grepl(x = colnames(.), pattern = blacklist)] %>%
                mutate_if(is.character, as.factor) %>%
                select_if(is.factor) %>%
                mutate_all(droplevels) %>%
                unique() %>%
                camel()
            assert_has_no_duplicates(data[["sampleName"]])
            rownames(data) <- data[["sampleID"]]
            data
        }
        if (is.character(interestingGroups)) {
            data <- uniteInterestingGroups(data, interestingGroups)
        }
        data <- sanitizeSampleData(data)
        assertHasRownames(data)
        as(data, return)
    }
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData<-",
    signature(
        object = "seurat",
        value = "ANY"
    ),
    getMethod("sampleData<-", "SingleCellExperiment")
)
