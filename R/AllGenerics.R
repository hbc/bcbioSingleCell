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
#' @inheritParams AllGenerics
#' @export
setGeneric("aggregateFeatures", function(object, ...) {
    standardGeneric("aggregateFeatures")
})



#' @rdname calculateMetrics
#' @inheritParams AllGenerics
#' @export
setGeneric("calculateMetrics", function(object, ...) {
    standardGeneric("calculateMetrics")
})



#' @rdname fetchTSNEData
#' @inheritParams AllGenerics
#' @export
setGeneric("fetchTSNEData", function(object, ...) {
    standardGeneric("fetchTSNEData")
})



#' @rdname fetchTSNEExpressionData
#' @inheritParams AllGenerics
#' @export
setGeneric("fetchTSNEExpressionData", function(object, ...) {
    standardGeneric("fetchTSNEExpressionData")
})



#' @rdname filterCells
#' @inheritParams AllGenerics
#' @export
setGeneric("filterCells", function(object, ...) {
    standardGeneric("filterCells")
})



#' @rdname knownMarkersDetected
#' @inheritParams AllGenerics
#' @export
setGeneric("knownMarkersDetected", function(object, knownMarkers, ...) {
    standardGeneric("knownMarkersDetected")
})



#' @rdname pcCutoff
#' @inheritParams AllGenerics
#' @export
setGeneric("pcCutoff", function(object, ...) {
    standardGeneric("pcCutoff")
})



#' @rdname plotCellCounts
#' @inheritParams AllGenerics
#' @export
setGeneric("plotCellCounts", function(object, ...) {
    standardGeneric("plotCellCounts")
})



#' @rdname plotFeatures
#' @inheritParams AllGenerics
#' @export
setGeneric("plotFeatures", function(object, ...) {
    standardGeneric("plotFeatures")
})



#' @rdname plotGenesPerCell
#' @inheritParams AllGenerics
#' @export
setGeneric("plotGenesPerCell", function(object, ...) {
    standardGeneric("plotGenesPerCell")
})



#' @rdname plotKnownMarkers
#' @inheritParams AllGenerics
#' @export
setGeneric("plotKnownMarkers", function(object, knownMarkers, ...) {
    standardGeneric("plotKnownMarkers")
})



#' @rdname plotMarkers
#' @inheritParams AllGenerics
#' @export
setGeneric("plotMarkers", function(object, ...) {
    standardGeneric("plotMarkers")
})



#' @rdname plotMitoRatio
#' @inheritParams AllGenerics
#' @export
setGeneric("plotMitoRatio", function(object, ...) {
    standardGeneric("plotMitoRatio")
})



#' @rdname plotNovelty
#' @inheritParams AllGenerics
#' @export
setGeneric("plotNovelty", function(object, ...) {
    standardGeneric("plotNovelty")
})



#' @rdname plotReadsPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#' @inheritParams AllGenerics
#' @export
setGeneric("plotReadsPerCell", function(object, ...) {
    standardGeneric("plotReadsPerCell")
})



#' @rdname plotStressGenes
#' @inheritParams AllGenerics
#' @export
setGeneric("plotStressGenes", function(object, ...) {
    standardGeneric("plotStressGenes")
})



#' @rdname plotTopMarkers
#' @inheritParams AllGenerics
#' @export
setGeneric("plotTopMarkers", function(object, topMarkers, ...) {
    standardGeneric("plotTopMarkers")
})



#' @rdname plotUMIsPerCell
#' @inheritParams AllGenerics
#' @export
setGeneric("plotUMIsPerCell", function(object, ...) {
    standardGeneric("plotUMIsPerCell")
})



#' @rdname plotUMIsVsGenes
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#' @inheritParams AllGenerics
#' @inherit plotGenesPerCell
#' @export
setGeneric("plotUMIsVsGenes", function(object, ...) {
    standardGeneric("plotUMIsVsGenes")
})



#' @rdname prepareSingleCellTemplate
#' @inheritParams AllGenerics
#' @export
setGeneric("prepareSingleCellTemplate", function(object, ...) {
    standardGeneric("prepareSingleCellTemplate")
})



#' @rdname plotZerosVsDepth
#' @inheritParams AllGenerics
#' @export
setGeneric("plotZerosVsDepth", function(object, ...) {
    standardGeneric("plotZerosVsDepth")
})



#' @rdname quantileHeatmap
#' @inheritParams AllGenerics
#' @export
setGeneric("quantileHeatmap", function(object, ...) {
    standardGeneric("quantileHeatmap")
})



#' @rdname readMarkers
#' @inheritParams AllGenerics
#' @export
setGeneric("readMarkers", function(object, ...) {
    standardGeneric("readMarkers")
})



#' @rdname subsetPerSample
#' @inheritParams AllGenerics
#' @export
setGeneric("subsetPerSample", function(object, ...) {
    standardGeneric("subsetPerSample")
})



#' @rdname topBarcodes
#' @inheritParams AllGenerics
#' @export
setGeneric("topBarcodes", function(object, ...) {
    standardGeneric("topBarcodes")
})



#' @rdname topMarkers
#' @inheritParams AllGenerics
#' @export
setGeneric("topMarkers", function(object, ...) {
    standardGeneric("topMarkers")
})
