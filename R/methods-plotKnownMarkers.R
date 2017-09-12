#' Plot Known Markers
#'
#' @rdname plotKnownMarkers
#' @name plotKnownMarkers
#' @inherit plotMarkers
#'
#' @param knownMarkers [knownMarkersDetected()] [tibble] grouped by cluster.
NULL



# Constructor ====
.plotKnownMarkers <- function(object, knownMarkers, headerLevel = 2L) {
    if (nrow(knownMarkers) == 0L) {
        return(NULL)
    }
    cellTypes <- knownMarkers %>%
        pull("cell") %>%
        unique
    pblapply(seq_along(cellTypes), function(a) {
        cellType <- cellTypes[[a]]
        symbols <- knownMarkers %>%
            tidy_filter(.data[["cell"]] == !!cellType) %>%
            pull("symbol") %>%
            unique %>%
            sort
        if (!is.null(symbols)) {
            message(cellType)
            mdHeader(cellType, level = headerLevel, tabset = TRUE, asis = TRUE)
            plotMarkers(object, symbols, headerLevel = headerLevel + 1L)
        } else {
            NULL
        }
    }) %>%
        invisible
}



# Methods ====
#' @rdname plotKnownMarkers
#' @export
setMethod(
    "plotKnownMarkers",
    signature(object = "seurat", knownMarkers = "grouped_df"),
    .plotKnownMarkers)
