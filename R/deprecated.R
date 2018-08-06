# nocov start



#' Deprecated Functions
#'
#' @name deprecated
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return [.Deprecated()].
NULL



#' Defunct Functions
#'
#' @name defunct
#' @keywords internal
#'
#' @inheritParams general
#'
#' @return [.Defunct()].
NULL



# v0.0.19 ======================================================================
#' @rdname defunct
#' @export
loadSingleCellRun <- function(...) {
    .Defunct("bcbioSingleCell")
}



# v0.0.24 ======================================================================
#' @rdname defunct
#' @export
darkTheme <- function(...) {
    .Defunct("basejump::theme_midnight")
}



#' @rdname defunct
#' @export
pcCutoff <- function(...) {
    .Defunct("plotPCElbow")
}



# v0.1.0 =======================================================================
#' @rdname deprecated
#' @export
calculateMetrics <- function(...) {
    .Deprecated("metrics")
    metrics(...)
}



# v0.1.1 =======================================================================
#' @rdname deprecated
#' @export
inflectionPoint <- function(...) {
    .Defunct("plotBarcodeRanks")
}



#' @rdname deprecated
#' @export
plotCumulativeUMIsPerCell <- function(...) {
    .Defunct("plotUMIsPerCell")
}



# v0.1.2 =======================================================================
#' @rdname deprecated
#' @export
loadCellRanger <- function(...) {
    .Deprecated("readCellRanger")
    readCellRanger(...)
}



#' @rdname deprecated
#' @export
loadSingleCell <- function(...) {
    .Deprecated("bcbioSingleCell")
    bcbioSingleCell(...)
}



# v0.2.0 =======================================================================
#' @rdname deprecated
#' @export
cellCountsPerCluster <- function(...) {
    .Deprecated("pointillism::cellCountsPerCluster")
    requireNamespace("pointillism")
    pointillism::cellCountsPerCluster(...)
}

#' @rdname deprecated
#' @export
cellTypesPerCluster <- function(...) {
    .Deprecated("pointillism::cellTypesPerCluster")
    requireNamespace("pointillism")
    pointillism::cellTypesPerCluster(...)
}

#' @rdname deprecated
#' @export
clusterCellCountsPerSample <- function(...) {
    .Deprecated("pointillism::clusterCellCountsPerSample")
    requireNamespace("pointillism")
    pointillism::clusterCellCountsPerSample(...)
}

#' @rdname deprecated
#' @export
diffExp <- function(...) {
    .Deprecated("pointillism::diffExp")
    requireNamespace("pointillism")
    pointillism::diffExp(...)
}

#' @rdname deprecated
#' @export
fetchGeneData <- function(...) {
    .Deprecated("pointillism::fetchGeneData")
    requireNamespace("pointillism")
    pointillism::fetchGeneData(...)
}

#' @rdname deprecated
#' @export
fetchPCAData <- function(...) {
    .Deprecated("pointillism::fetchPCAData")
    requireNamespace("pointillism")
    pointillism::fetchPCAData(...)
}

#' @rdname deprecated
#' @export
fetchTSNEData <- function(...) {
    .Deprecated("pointillism::fetchTSNEData")
    requireNamespace("pointillism")
    pointillism::fetchTSNEData(...)
}

#' @rdname deprecated
#' @export
fetchTSNEExpressionData <- function(...) {
    .Deprecated("pointillism::fetchTSNEExpressionData")
    requireNamespace("pointillism")
    pointillism::fetchTSNEExpressionData(...)
}

#' @rdname deprecated
#' @export
fetchUMAPData <- function(...) {
    .Deprecated("pointillism::fetchUMAPData")
    requireNamespace("pointillism")
    pointillism::fetchUMAPData(...)
}

#' @rdname deprecated
#' @export
fetchUMAPExpressionData <- function(...) {
    .Deprecated("pointillism::fetchUMAPExpressionData")
    requireNamespace("pointillism")
    pointillism::fetchUMAPExpressionData(...)
}

#' @rdname deprecated
#' @export
knownMarkersDetected <- function(...) {
    .Deprecated("pointillism::knownMarkersDetected")
    requireNamespace("pointillism")
    pointillism::knownMarkersDetected(...)
}

#' @rdname deprecated
#' @export
plotCellTypesPerCluster <- function(...) {
    .Deprecated("pointillism::plotCellTypesPerCluster")
    requireNamespace("pointillism")
    pointillism::plotCellTypesPerCluster(...)
}

#' @rdname deprecated
#' @export
plotFeatureTSNE <- function(...) {
    .Deprecated("pointillism::plotFeatureTSNE")
    requireNamespace("pointillism")
    pointillism::plotFeatureTSNE(...)
}

#' @rdname deprecated
#' @export
plotFeatureUMAP <- function(...) {
    .Deprecated("pointillism::plotFeatureUMAP")
    requireNamespace("pointillism")
    pointillism::plotFeatureUMAP(...)
}

#' @rdname deprecated
#' @export
plotGene <- function(...) {
    .Deprecated("pointillism::plotGene")
    requireNamespace("pointillism")
    pointillism::plotGene(...)
}

#' @rdname deprecated
#' @export
plotKnownMarkers <- function(...) {
    .Deprecated("pointillism::plotKnownMarkers")
    requireNamespace("pointillism")
    pointillism::plotKnownMarkers(...)
}

#' @rdname deprecated
#' @export
plotKnownMarkersDetected <- function(...) {
    .Deprecated("pointillism::plotKnownMarkersDetected")
    requireNamespace("pointillism")
    pointillism::plotKnownMarkersDetected(...)
}

#' @rdname deprecated
#' @export
plotMarker <- function(...) {
    .Deprecated("pointillism::plotMarker")
    requireNamespace("pointillism")
    pointillism::plotMarker(...)
}

#' @rdname deprecated
#' @export
plotMarkerTSNE <- function(...) {
    .Deprecated("pointillism::plotMarkerTSNE")
    requireNamespace("pointillism")
    pointillism::plotMarkerTSNE(...)
}

#' @rdname deprecated
#' @export
plotMarkerUMAP <- function(...) {
    .Deprecated("pointillism::plotMarkerUMAP")
    requireNamespace("pointillism")
    pointillism::plotMarkerUMAP(...)
}

#' @rdname deprecated
#' @export
plotMarkers <- function(...) {
    .Deprecated("pointillism::plotMarkers")
    requireNamespace("pointillism")
    pointillism::plotMarkers(...)
}

#' @rdname deprecated
#' @export
plotPCElbow <- function(...) {
    .Deprecated("pointillism::plotPCElbow")
    requireNamespace("pointillism")
    pointillism::plotPCElbow(...)
}

#' @rdname deprecated
#' @export
plotTSNE <- function(...) {
    .Deprecated("pointillism::plotTSNE")
    requireNamespace("pointillism")
    pointillism::plotTSNE(...)
}

#' @rdname deprecated
#' @export
plotTSNEExpressionData <- function(...) {
    .Deprecated("pointillism::plotTSNEExpressionData")
    requireNamespace("pointillism")
    pointillism::plotTSNEExpressionData(...)
}

#' @rdname deprecated
#' @export
plotTopMarkers <- function(...) {
    .Deprecated("pointillism::plotTopMarkers")
    requireNamespace("pointillism")
    pointillism::plotTopMarkers(...)
}

#' @rdname deprecated
#' @export
plotUMAP <- function(...) {
    .Deprecated("pointillism::plotUMAP")
    requireNamespace("pointillism")
    pointillism::plotUMAP(...)
}

#' @rdname deprecated
#' @export
readCellTypeMarkers <- function(...) {
    .Deprecated("pointillism::readCellTypeMarkers")
    requireNamespace("pointillism")
    pointillism::readCellTypeMarkers(...)
}

#' @rdname deprecated
#' @export
sanitizeMarkers <- function(...) {
    .Deprecated("pointillism::sanitizeSeuratMarkers")
    requireNamespace("pointillism")
    pointillism::sanitizeSeuratMarkers(...)
}

#' @rdname deprecated
#' @export
topMarkers <- function(...) {
    .Deprecated("pointillism::topMarkers")
    requireNamespace("pointillism")
    pointillism::topMarkers(...)
}



# nocov end
