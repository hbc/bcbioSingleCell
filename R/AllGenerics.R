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
#' @param ... Additional arguments.
#'
#' @return No value.
NULL



#' @rdname aggregateFeatures
#' @family Data Management Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("aggregateFeatures", function(object, ...) {
    standardGeneric("aggregateFeatures")
})



#' @rdname aggregateReplicates
#' @family Data Management Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("aggregateReplicates", function(object, ...) {
    standardGeneric("aggregateReplicates")
})



#' @rdname bcbio
#' @family `bcbio` Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("bcbio", function(object, ...) {
    standardGeneric("bcbio")
})



#' @rdname bcbio
#' @inheritParams AllGenerics
#' @export
setGeneric("bcbio<-", function(object, ..., value) {
    standardGeneric("bcbio<-")
})



#' @rdname calculateMetrics
#' @family QC Metrics Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("calculateMetrics", function(object, ...) {
    standardGeneric("calculateMetrics")
})



#' @rdname downloads
#' @family `bcbio` Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("downloads", function(object) {
    standardGeneric("downloads")
})



#' @rdname filterCells
#' @family Data Management Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("filterCells", function(object, ...) {
    standardGeneric("filterCells")
})



#' @rdname loadSinglecell
#' @family `bcbio` Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("loadSinglecell", function(object, ...) {
    standardGeneric("loadSinglecell")
})



#' @rdname metadataTable
#' @family RMarkdown Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("metadataTable", function(object, ...) {
    standardGeneric("metadataTable")
})



#' @rdname interestingGroups
#' @family Data Management Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("interestingGroups", function(object) {
    standardGeneric("interestingGroups")
})



#' @rdname knownMarkersDetected
#' @family Clustering Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("knownMarkersDetected", function(object, known, ...) {
    standardGeneric("knownMarkersDetected")
})



#' @rdname metrics
#' @family Quality Control Metrics
#' @inheritParams AllGenerics
#' @export
setGeneric("metrics", function(object, ...) {
    standardGeneric("metrics")
})



#' @rdname plotCellCounts
#' @family Quality Control Metrics
#' @inheritParams AllGenerics
#' @export
setGeneric("plotCellCounts", function(object, ...) {
    standardGeneric("plotCellCounts")
})



#' @rdname plotCellularBarcodes
#' @family Quality Control Metrics
#' @inheritParams AllGenerics
#' @export
setGeneric("plotCellularBarcodes", function(object, ...) {
    standardGeneric("plotCellularBarcodes")
})



#' @rdname plotClusters
#' @family Clustering Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("plotClusters", function(object, ...) {
    standardGeneric("plotClusters")
})



#' @rdname plotGenesPerCell
#' @family Quality Control Metrics
#' @inheritParams AllGenerics
#' @export
setGeneric("plotGenesPerCell", function(object, ...) {
    standardGeneric("plotGenesPerCell")
})



#' @rdname plotKnownMarkers
#' @family Clustering Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("plotKnownMarkers", function(x, y, ...) {
    standardGeneric("plotKnownMarkers")
})



#' @rdname plotMitoRatio
#' @family Quality Control Metrics
#' @inheritParams AllGenerics
#' @export
setGeneric("plotMitoRatio", function(object, ...) {
    standardGeneric("plotMitoRatio")
})



#' @rdname plotNovelty
#' @family Quality Control Metrics
#' @inheritParams AllGenerics
#' @export
setGeneric("plotNovelty", function(object, ...) {
    standardGeneric("plotNovelty")
})



#' @rdname plotStressGenes
#' @family Clustering Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("plotStressGenes", function(object, ...) {
    standardGeneric("plotStressGenes")
})



#' @rdname plotTopMarkers
#' @family Clustering Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("plotTopMarkers", function(object, markers, ...) {
    standardGeneric("plotTopMarkers")
})



#' @rdname plotUMIsPerCell
#' @family Quality Control Metrics
#' @inheritParams AllGenerics
#' @export
setGeneric("plotUMIsPerCell", function(object, ...) {
    standardGeneric("plotUMIsPerCell")
})



#' @rdname plotUMIsVsGenes
#' @family Quality Control Metrics
#' @inheritParams AllGenerics
#' @export
setGeneric("plotUMIsVsGenes", function(object, ...) {
    standardGeneric("plotUMIsVsGenes")
})



#' @rdname readKnownMarkers
#' @family Data Management Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("readKnownMarkers", function(object, ...) {
    standardGeneric("readKnownMarkers")
})



#' @rdname sampleMetadata
#' @family Data Management Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("sampleMetadata", function(object, ...) {
    standardGeneric("sampleMetadata")
})



#' @rdname selectSamples
#' @family Data Management Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("selectSamples", function(object, ...) {
    standardGeneric("selectSamples")
})



#' @rdname topBarcodes
#' @family Clustering Utilities
#' @inheritParams AllGenerics
#' @export
setGeneric("topBarcodes", function(object, ...) {
    standardGeneric("topBarcodes")
})



#' @rdname topMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("topMarkers", function(object, ...) {
    standardGeneric("topMarkers")
})
