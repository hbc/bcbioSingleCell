#' Plot Gene Markers
#'
#' @rdname plotMarkers
#' @name plotMarkers
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
        .plotFeatures.seurat(object, symbols, nCol = 2L)
    } else {
        lapply(seq_along(symbols), function(a) {
            mdHeader(symbols[[a]], level = headerLevel, asis = TRUE)
            .plotFeatures.seurat(object, symbols[[a]])
        }) %>%
            invisible
    }
})
