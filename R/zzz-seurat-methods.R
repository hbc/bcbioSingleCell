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
    function(x, ...) {
        assay(.as.SingleCellExperiment.seurat(x), ...)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assayNames
#' @export
setMethod(
    "assayNames",
    signature("seurat"),
    function(x, ...) {
        assayNames(.as.SingleCellExperiment.seurat(x), ...)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment assays
#' @export
setMethod(
    "assays",
    signature("seurat"),
    function(x, ...) {
        assays(.as.SingleCellExperiment.seurat(x), ...)
    }
)



#' @rdname barcodeRanksPerSample
#' @export
setMethod(
    "barcodeRanksPerSample",
    signature("seurat"),
    getMethod("barcodeRanksPerSample", "SingleCellExperiment")
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "cell2sample",
    signature("seurat"),
    getMethod("cell2sample", "SingleCellExperiment")
)



#' @rdname cellCountsPerCluster
#' @export
setMethod(
    "cellCountsPerCluster",
    signature("seurat"),
    getMethod("cellCountsPerCluster", "SingleCellExperiment")
)



#' @rdname clusterCellCountsPerSample
#' @export
setMethod(
    "clusterCellCountsPerSample",
    signature("seurat"),
    getMethod("clusterCellCountsPerSample", "SingleCellExperiment")
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment colData
#' @export
setMethod(
    "colData",
    signature("seurat"),
    function(x, ...) {
        colData(.as.SingleCellExperiment.seurat(x), ...)
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
        colnames(.as.SingleCellExperiment.seurat(x))
    }
)



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
#' @export
setMethod(
    "counts",
    signature("seurat"),
    function(object, ...) {
        counts(.as.SingleCellExperiment.seurat(object), ...)
    }
)



#' @rdname diffExp
#' @export
setMethod(
    "diffExp",
    signature("seurat"),
    getMethod("diffExp", "SingleCellExperiment")
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchGeneData",
    signature("seurat"),
    getMethod("fetchGeneData", "SingleCellExperiment")
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchReducedDimData",
    signature("seurat"),
    getMethod("fetchReducedDimData", "SingleCellExperiment")
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchReducedDimExpressionData",
    signature("seurat"),
    getMethod("fetchReducedDimExpressionData", "SingleCellExperiment")
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
#' @importFrom S4Vectors metadata
#' @export
setMethod(
    "metadata",
    signature("seurat"),
    function(x, ...) {
        stash <- slot(x, "misc")[["bcbio"]][["metadata"]]
        if (!is.null(stash)) {
            return(stash)
        }
        metadata(.as.SingleCellExperiment.seurat(x), ...)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom S4Vectors metadata<-
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



#' @rdname metrics
#' @export
setMethod(
    "metrics",
    signature("seurat"),
    function(object, ...) {
        metrics(as(object, "SingleCellExperiment"))
    }
)



#' @rdname metricsPerSample
#' @export
setMethod(
    "metricsPerSample",
    signature("seurat"),
    getMethod("metricsPerSample", "SingleCellExperiment")
)



#' @rdname plotBarcodeRanks
#' @export
setMethod(
    "plotBarcodeRanks",
    signature("seurat"),
    getMethod("plotBarcodeRanks", "SingleCellExperiment")
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "plotCellCounts",
    signature("seurat"),
    getMethod("plotCellCounts", "SingleCellExperiment")
)



#' @rdname plotCellTypesPerCluster
#' @export
setMethod(
    "plotCellTypesPerCluster",
    signature("seurat"),
    getMethod("plotCellTypesPerCluster", "SingleCellExperiment")
)



#' @rdname plotDot
#' @export
setMethod(
    "plotDot",
    signature("seurat"),
    getMethod("plotDot", "SingleCellExperiment")
)



#' @rdname plotFeature
#' @export
setMethod(
    "plotFeature",
    signature("seurat"),
    getMethod("plotFeature", "SingleCellExperiment")
)



#' @rdname plotGene
#' @export
setMethod(
    "plotGene",
    signature("seurat"),
    getMethod("plotGene", "SingleCellExperiment")
)



#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("seurat"),
    getMethod("plotGenesPerCell", "SingleCellExperiment")
)



#' @rdname plotMarker
#' @export
setMethod(
    "plotKnownMarkersDetected",
    signature("seurat"),
    getMethod("plotKnownMarkersDetected", "SingleCellExperiment")
)



#' @rdname plotMarker
#' @export
setMethod(
    "plotMarker",
    signature("seurat"),
    getMethod("plotMarker", "SingleCellExperiment")
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



#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    signature("seurat"),
    getMethod("plotNovelty", "SingleCellExperiment")
)



#' @rdname plotReducedDim
#' @export
setMethod(
    "plotPCA",
    signature("seurat"),
    getMethod("plotPCA", "SingleCellExperiment")
)



#' @rdname plotQC
#' @export
setMethod(
    "plotQC",
    signature("seurat"),
    getMethod("plotQC", "SingleCellExperiment")
)



#' @rdname plotReducedDim
#' @export
setMethod(
    "plotReducedDim",
    signature("seurat"),
    getMethod("plotReducedDim", "SingleCellExperiment")
)



#' @rdname plotMarker
#' @export
setMethod(
    "plotTopMarkers",
    signature("seurat"),
    getMethod("plotTopMarkers", "SingleCellExperiment")
)



#' @rdname plotReducedDim
#' @export
setMethod(
    "plotTSNE",
    signature("seurat"),
    getMethod("plotTSNE", "SingleCellExperiment")
)



#' @rdname plotReducedDim
#' @export
setMethod(
    "plotUMAP",
    signature("seurat"),
    getMethod("plotUMAP", "SingleCellExperiment")
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



#' @rdname plotViolin
#' @export
setMethod(
    "plotViolin",
    signature("seurat"),
    getMethod("plotViolin", "SingleCellExperiment")
)



#' @rdname plotZerosVsDepth
#' @export
setMethod(
    "plotZerosVsDepth",
    signature("seurat"),
    getMethod("plotZerosVsDepth", "SingleCellExperiment")
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SingleCellExperiment reducedDims
#' @export
setMethod(
    "reducedDims",
    signature("seurat"),
    function(x) {
        reducedDims(.as.SingleCellExperiment.seurat(x))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment rowData
#' @export
setMethod(
    "rowData",
    signature("seurat"),
    function(x, ...) {
        rowData(as(x, "SingleCellExperiment"), ...)
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom BiocGenerics rownames
#' @export
setMethod(
    "rownames",
    signature("seurat"),
    function(x) {
        rownames(.as.SingleCellExperiment.seurat(x))
    }
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom SummarizedExperiment rowRanges
#' @export
setMethod(
    "rowRanges",
    signature("seurat"),
    function(x, ...) {
        gr <- rowRanges(.as.SingleCellExperiment.seurat(x), ...)

        # Attempt to use stashed rowRanges, if present
        stash <- slot(x, "misc")[["bcbio"]][["rowRanges"]]
        if (is(stash, "GRanges")) {
            assert_is_subset(c("geneID", "geneName"), colnames(mcols(stash)))
            # Check to see if we're using IDs or symbols
            if (any(names(gr) %in% mcols(stash)[["geneID"]])) {
                col <- "geneID"
            } else if (any(names(gr) %in% mcols(stash)[["geneName"]])) {
                col <- "geneName"
            } else {
                stop("Failed to match identifiers to rownames")
            }
            names(stash) <- make.unique(as.character(mcols(stash)[[col]]))
            assert_is_subset(names(gr), names(stash))
            stash <- stash[names(gr)]
            assert_are_disjoint_sets(
                x = colnames(mcols(gr)),
                y = colnames(mcols(stash))
            )
            mcols(stash) <- cbind(mcols(stash), mcols(gr))
            gr <- stash
        }

        gr
    }
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData",
    signature("seurat"),
    getMethod("sampleData", "SingleCellExperiment")
)



#' @rdname sampleData
#' @export
setMethod(
    "sampleData<-",
    signature(
        object = "seurat",
        value = "DataFrame"
    ),
    getMethod(
        "sampleData<-",
        signature(
            object = "SingleCellExperiment",
            value = "DataFrame"
        )
    )
)



#' @rdname seurat-SingleCellExperiment
#' @importFrom basejump sampleNames
#' @export
setMethod(
    "sampleNames",
    signature("seurat"),
    getMethod("sampleNames", "SummarizedExperiment")
)



#' @rdname topBarcodes
#' @export
setMethod(
    "topBarcodes",
    signature("seurat"),
    getMethod("topBarcodes", "SingleCellExperiment")
)
