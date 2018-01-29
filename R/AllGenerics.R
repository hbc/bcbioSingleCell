#' S4 Generics
#'
#' @rdname AllGenerics
#' @name AllGenerics
#' @keywords internal
#'
#' @param object Object.
#' @param x Primary object.
#' @param y Secondary object.
#' @param value Object to assign.
#' @param ... *Additional arguments (for the S4 generic definition).*
#'
#' @return No value.
NULL



#' @rdname calculateMetrics
#' @export
setGeneric("calculateMetrics", function(object, ...) {
    standardGeneric("calculateMetrics")
})



#' @rdname cell2sample
#' @export
setGeneric("cell2sample", function(object, ...) {
    standardGeneric("cell2sample")
})



#' @rdname cellTypesPerCluster
#' @export
setGeneric("cellTypesPerCluster", function(object, ...) {
    standardGeneric("cellTypesPerCluster")
})



#' @rdname fetchGeneData
#' @export
setGeneric("fetchGeneData", function(object, ...) {
    standardGeneric("fetchGeneData")
})



#' @rdname fetchPCAData
#' @export
setGeneric("fetchPCAData", function(object, ...) {
    standardGeneric("fetchPCAData")
})



#' @rdname fetchTSNEData
#' @export
setGeneric("fetchTSNEData", function(object, ...) {
    standardGeneric("fetchTSNEData")
})



#' @rdname fetchTSNEExpressionData
#' @export
setGeneric("fetchTSNEExpressionData", function(object, ...) {
    standardGeneric("fetchTSNEExpressionData")
})



#' @rdname filterCells
#' @export
setGeneric("filterCells", function(object, ...) {
    standardGeneric("filterCells")
})



#' @rdname knownMarkersDetected
#' @export
setGeneric("knownMarkersDetected", function(all, known, ...) {
    standardGeneric("knownMarkersDetected")
})



#' @rdname plotCellCounts
#' @export
setGeneric("plotCellCounts", function(object, ...) {
    standardGeneric("plotCellCounts")
})



#' @rdname plotCellTypesPerCluster
#' @export
setGeneric(
    "plotCellTypesPerCluster",
    function(object, cellTypesPerCluster, ...) {
        standardGeneric("plotCellTypesPerCluster")
    })



#' @rdname plotFeatureTSNE
#' @export
setGeneric("plotFeatureTSNE", function(object, ...) {
    standardGeneric("plotFeatureTSNE")
})



#' @rdname plotGenesPerCell
#' @export
setGeneric("plotGenesPerCell", function(object, ...) {
    standardGeneric("plotGenesPerCell")
})



#' @rdname plotKnownMarkersDetected
#' @export
setGeneric(
    "plotKnownMarkersDetected",
    function(object, knownMarkersDetected, ...) {
        standardGeneric("plotKnownMarkersDetected")
    })



#' @rdname plotMarkers
#' @export
setGeneric("plotMarkers", function(object, ...) {
    standardGeneric("plotMarkers")
})



#' @rdname plotMarkerTSNE
#' @export
setGeneric("plotMarkerTSNE", function(object, ...) {
    standardGeneric("plotMarkerTSNE")
})



#' @rdname plotMitoRatio
#' @export
setGeneric("plotMitoRatio", function(object, ...) {
    standardGeneric("plotMitoRatio")
})



#' @rdname plotMitoVsCoding
#' @export
setGeneric("plotMitoVsCoding", function(object, ...) {
    standardGeneric("plotMitoVsCoding")
})



#' @rdname plotNovelty
#' @export
setGeneric("plotNovelty", function(object, ...) {
    standardGeneric("plotNovelty")
})



#' @rdname plotPCElbow
#' @export
setGeneric("plotPCElbow", function(object, ...) {
    standardGeneric("plotPCElbow")
})



#' @rdname plotReadsPerCell
#' @export
setGeneric("plotReadsPerCell", function(object, ...) {
    standardGeneric("plotReadsPerCell")
})



#' @rdname plotTopMarkers
#' @export
setGeneric("plotTopMarkers", function(object, topMarkers, ...) {
    standardGeneric("plotTopMarkers")
})



#' @rdname plotTSNE
#' @export
setGeneric("plotTSNE", function(object, ...) {
    standardGeneric("plotTSNE")
})



#' @rdname plotUMIsPerCell
#' @export
setGeneric("plotUMIsPerCell", function(object, ...) {
    standardGeneric("plotUMIsPerCell")
})



#' @rdname plotUMIsVsGenes
#' @export
setGeneric("plotUMIsVsGenes", function(object, ...) {
    standardGeneric("plotUMIsVsGenes")
})



#' @rdname plotZerosVsDepth
#' @export
setGeneric("plotZerosVsDepth", function(object, ...) {
    standardGeneric("plotZerosVsDepth")
})



#' @rdname prepareSingleCellTemplate
#' @export
setGeneric("prepareSingleCellTemplate", function(object, ...) {
    standardGeneric("prepareSingleCellTemplate")
})



#' @rdname readCellTypeMarkersFile
#' @export
setGeneric("readCellTypeMarkersFile", function(object, ...) {
    standardGeneric("readCellTypeMarkersFile")
})



#' @rdname sanitizeMarkers
#' @export
setGeneric("sanitizeMarkers", function(object, markers, ...) {
    standardGeneric("sanitizeMarkers")
})



#' @rdname subsetPerSample
#' @export
setGeneric("subsetPerSample", function(object, ...) {
    standardGeneric("subsetPerSample")
})



#' @rdname topBarcodes
#' @export
setGeneric("topBarcodes", function(object, ...) {
    standardGeneric("topBarcodes")
})



#' @rdname topMarkers
#' @export
setGeneric("topMarkers", function(object, ...) {
    standardGeneric("topMarkers")
})
