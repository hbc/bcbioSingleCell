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
#' @author Rory Kirchner, Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("aggregateFeatures", function(object, ...) {
    standardGeneric("aggregateFeatures")
})



#' @rdname aggregateReplicates
#' @family Data Management Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("aggregateReplicates", function(object, ...) {
    standardGeneric("aggregateReplicates")
})



#' @rdname bcbio
#' @author Michael Steinbaugh
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
#' @author Rory Kirchner, Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("calculateMetrics", function(object, ...) {
    standardGeneric("calculateMetrics")
})



#' @rdname download
#' @family `bcbio` Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("download", function(object) {
    standardGeneric("download")
})



#' @rdname filterCells
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("filterCells", function(object, ...) {
    standardGeneric("filterCells")
})



#' @rdname interestingGroups
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("interestingGroups", function(object) {
    standardGeneric("interestingGroups")
})



#' @rdname knownMarkersDetected
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("knownMarkersDetected", function(object, knownMarkers, ...) {
    standardGeneric("knownMarkersDetected")
})



#' @rdname loadSingleCellRun
#' @author Michael Steinbaugh, Rory Kirchner
#' @inheritParams AllGenerics
#' @export
setGeneric("loadSingleCellRun", function(object, ...) {
    standardGeneric("loadSingleCellRun")
})



#' @rdname metadataTable
#' @family RMarkdown Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("metadataTable", function(object, ...) {
    standardGeneric("metadataTable")
})



#' @rdname metrics
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#' @inheritParams AllGenerics
#' @export
setGeneric("metrics", function(object, ...) {
    standardGeneric("metrics")
})



#' @rdname plotCellCounts
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#' @inheritParams AllGenerics
#' @inherit plotGenesPerCell
#' @export
setGeneric("plotCellCounts", function(object, ...) {
    standardGeneric("plotCellCounts")
})



#' @rdname plotClusters
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotClusters", function(object, ...) {
    standardGeneric("plotClusters")
})



#' @rdname plotGenesPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#' @inheritParams AllGenerics
#' @export
setGeneric("plotGenesPerCell", function(object, ...) {
    standardGeneric("plotGenesPerCell")
})



#' @rdname plotKnownMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotKnownMarkers", function(object, knownMarkers, ...) {
    standardGeneric("plotKnownMarkers")
})



#' @rdname plotMitoRatio
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#' @inheritParams AllGenerics
#' @inherit plotGenesPerCell
#' @export
setGeneric("plotMitoRatio", function(object, ...) {
    standardGeneric("plotMitoRatio")
})



#' @rdname plotNovelty
#' @family Quality Control Metrics
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @inherit plotGenesPerCell
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
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotStressGenes", function(object, ...) {
    standardGeneric("plotStressGenes")
})



#' @rdname plotTopMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotTopMarkers", function(object, topMarkers, ...) {
    standardGeneric("plotTopMarkers")
})



#' @rdname plotUMIsPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#' @inheritParams AllGenerics
#' @inherit plotGenesPerCell
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



#' @rdname plotZerosVsDepth
#' @family Quality Control Metrics
#' @author Rory Kirchner, Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("plotZerosVsDepth", function(object, ...) {
    standardGeneric("plotZerosVsDepth")
})



#' @rdname readMarkers
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("readMarkers", function(object, ...) {
    standardGeneric("readMarkers")
})



#' @rdname sampleMetadata
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("sampleMetadata", function(object, ...) {
    standardGeneric("sampleMetadata")
})



#' @rdname selectSamples
#' @family Data Management Utilities
#' @author Michael Steinbaugh
#' @inheritParams AllGenerics
#' @export
setGeneric("selectSamples", function(object, ...) {
    standardGeneric("selectSamples")
})



#' @rdname topBarcodes
#' @family Clustering Utilities
#' @author Michael Steinbaugh
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
