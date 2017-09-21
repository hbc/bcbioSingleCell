#' Plot Known Markers
#'
#' @rdname plotKnownMarkers
#' @name plotKnownMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inherit plotMarkers
#'
#' @param knownMarkers [knownMarkersDetected()] [tibble] grouped by cluster.
NULL



# Constructors ====
.plotKnownMarkers <- function(
    object,
    knownMarkers,
    headerLevel = 2L) {
    if (nrow(knownMarkers) == 0L) {
        return(NULL)
    }
    cellTypes <- knownMarkers %>%
        pull("cell") %>%
        unique
    pblapply(seq_along(cellTypes), function(a) {
        cellType <- cellTypes[[a]]
        symbols <- knownMarkers %>%
            dplyr::filter(.data[["cell"]] == !!cellType) %>%
            pull("symbol") %>%
            unique %>%
            sort
        if (!is.null(symbols)) {
            mdHeader(cellType, level = headerLevel, tabset = TRUE, asis = TRUE)
            plotMarkers(object,
                        symbols = symbols,
                        headerLevel = headerLevel + 1L)
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
