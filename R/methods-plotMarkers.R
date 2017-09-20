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
    headerLevel = 2L) {
    lapply(seq_along(symbols), function(a) {
        mdHeader(symbols[[a]], level = headerLevel, asis = TRUE)
        .plotFeatureSeurat(object, symbols[[a]])
    }) %>%
        invisible
})
