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



#' @rdname aggregateFeatures
#' @export
setGeneric("aggregateFeatures", function(object, ...) {
    standardGeneric("aggregateFeatures")
})



#' @rdname calculateMetrics
#' @export
setGeneric("calculateMetrics", function(object, ...) {
    standardGeneric("calculateMetrics")
})



#' @rdname cellTypesPerCluster
#' @export
setGeneric("cellTypesPerCluster", function(object, ...) {
    standardGeneric("cellTypesPerCluster")
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



#' @rdname pcCutoff
#' @export
setGeneric("pcCutoff", function(object, ...) {
    standardGeneric("pcCutoff")
})



#' @rdname plotCellCounts
#' @export
setGeneric("plotCellCounts", function(object, ...) {
    standardGeneric("plotCellCounts")
})



#' @rdname plotCellTypesPerCluster
#' @export
setGeneric("plotCellTypesPerCluster", function(object, ...) {
    standardGeneric("plotCellTypesPerCluster")
})



#' @rdname plotDot
#' @export
setGeneric("plotDot", function(object, ...) {
    standardGeneric("plotDot")
})



#' @rdname plotFeatures
#' @export
setGeneric("plotFeatures", function(object, ...) {
    standardGeneric("plotFeatures")
})



#' @rdname plotGenesPerCell
#' @export
setGeneric("plotGenesPerCell", function(object, ...) {
    standardGeneric("plotGenesPerCell")
})



#' @rdname plotKnownMarkers
#' @export
setGeneric("plotKnownMarkers", function(object, knownMarkers, ...) {
    standardGeneric("plotKnownMarkers")
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



#' @rdname plotNovelty
#' @export
setGeneric("plotNovelty", function(object, ...) {
    standardGeneric("plotNovelty")
})



#' @rdname plotReadsPerCell
#' @export
setGeneric("plotReadsPerCell", function(object, ...) {
    standardGeneric("plotReadsPerCell")
})



#' @rdname plotStressGenes
#' @export
setGeneric("plotStressGenes", function(object, ...) {
    standardGeneric("plotStressGenes")
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



#' @rdname prepareSingleCellTemplate
#' @export
setGeneric("prepareSingleCellTemplate", function(object, ...) {
    standardGeneric("prepareSingleCellTemplate")
})



#' @rdname plotZerosVsDepth
#' @export
setGeneric("plotZerosVsDepth", function(object, ...) {
    standardGeneric("plotZerosVsDepth")
})



#' @rdname quantileHeatmap
#' @export
setGeneric("quantileHeatmap", function(object, ...) {
    standardGeneric("quantileHeatmap")
})



#' @rdname readMarkers
#' @export
setGeneric("readMarkers", function(object, ...) {
    standardGeneric("readMarkers")
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
