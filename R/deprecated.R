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
#' @rdname defunct
#' @export
calculateMetrics <- function(...) {
    .Defunct("metrics")
}



# v0.1.1 =======================================================================
#' @rdname defunct
#' @export
inflectionPoint <- function(...) {
    .Defunct("plotBarcodeRanks")
}



#' @rdname defunct
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
#' @rdname defunct
#' @export
cellCountsPerCluster <- function(...) {
    .Defunct("pointillism::cellCountsPerCluster")
}

#' @rdname defunct
#' @export
cellTypesPerCluster <- function(...) {
    .Defunct("pointillism::cellTypesPerCluster")
}

#' @rdname defunct
#' @export
clusterCellCountsPerSample <- function(...) {
    .Defunct("pointillism::clusterCellCountsPerSample")
}

#' @rdname defunct
#' @export
diffExp <- function(...) {
    .Defunct("pointillism::diffExp")
}

#' @rdname defunct
#' @export
fetchGeneData <- function(...) {
    .Defunct()
}

#' @rdname defunct
#' @export
fetchPCAData <- function(...) {
    .Defunct()
}

#' @rdname defunct
#' @export
fetchTSNEData <- function(...) {
    .Defunct()
}

#' @rdname defunct
#' @export
fetchTSNEExpressionData <- function(...) {
    .Defunct()
}

#' @rdname defunct
#' @export
fetchUMAPData <- function(...) {
    .Defunct()
}

#' @rdname defunct
#' @export
fetchUMAPExpressionData <- function(...) {
    .Defunct()
}

#' @rdname defunct
#' @export
knownMarkersDetected <- function(...) {
    .Defunct("pointillism::knownMarkersDetected")
}

#' @rdname defunct
#' @export
plotCellTypesPerCluster <- function(...) {
    .Defunct("pointillism::plotCellTypesPerCluster")
}

#' @rdname defunct
#' @export
plotFeatureTSNE <- function(...) {
    .Defunct("pointillism::plotFeature")
}

#' @rdname defunct
#' @export
plotFeatureUMAP <- function(...) {
    .Defunct("pointillism::plotFeature")
}

#' @rdname defunct
#' @export
plotGene <- function(...) {
    .Defunct("pointillism::plotGene")
}

#' @rdname defunct
#' @export
plotKnownMarkers <- function(...) {
    .Defunct("pointillism::plotKnownMarkersDetected")
}

#' @rdname defunct
#' @export
plotKnownMarkersDetected <- function(...) {
    .Defunct("pointillism::plotKnownMarkersDetected")
}

#' @rdname defunct
#' @export
plotMarker <- function(...) {
    .Defunct("pointillism::plotMarker")
}

#' @rdname defunct
#' @export
plotMarkerTSNE <- function(...) {
    .Defunct("pointillism::plotMarker")
}

#' @rdname defunct
#' @export
plotMarkerUMAP <- function(...) {
    .Defunct("pointillism::plotMarker")
}

#' @rdname defunct
#' @export
plotMarkers <- function(...) {
    .Defunct("pointillism::plotMarkers")
}

#' @rdname defunct
#' @export
plotPCElbow <- function(...) {
    .Defunct("pointillism::plotPCElbow")
}

#' @rdname defunct
#' @export
plotTSNE <- function(...) {
    .Defunct("pointillism::plotTSNE")
}

#' @rdname defunct
#' @export
plotTSNEExpressionData <- function(...) {
    .Defunct("pointillism::plotMarker")
}

#' @rdname defunct
#' @export
plotTopMarkers <- function(...) {
    .Defunct("pointillism::plotTopMarkers")
}

#' @rdname defunct
#' @export
plotUMAP <- function(...) {
    .Defunct("pointillism::plotUMAP")
}

#' @rdname defunct
#' @export
readCellTypeMarkers <- function(...) {
    .Defunct("pointillism::readCellTypeMarkers")
}

#' @rdname defunct
#' @export
sanitizeMarkers <- function(...) {
    .Defunct("pointillism::sanitizeSeuratMarkers")
}

#' @rdname defunct
#' @export
topMarkers <- function(...) {
    .Defunct("pointillism::topMarkers")
}



# nocov end
