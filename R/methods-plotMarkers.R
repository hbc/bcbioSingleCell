#' Plot Gene Markers
#'
#' @rdname plotMarkers
#' @name plotMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inherit plotFeatures
#'
#' @param symbols Character vector of gene marker symbols.
NULL



# Methods ====
#' @rdname plotMarkers
#' @export
setMethod("plotMarkers", "seurat", function(
    object,
    symbols,
    headerLevel = 2L,
    combine = TRUE) {
    if (isTRUE(combine)) {
        .plotFeaturesSeurat(object, symbols, nCol = 2L)
    } else {
        lapply(seq_along(symbols), function(a) {
            mdHeader(symbols[[a]], level = headerLevel, asis = TRUE)
            .plotFeaturesSeurat(object, symbols[[a]])
        }) %>%
            invisible
    }
})
