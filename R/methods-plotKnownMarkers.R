#' Plot Known Markers
#'
#' @rdname plotKnownMarkers
#' @name plotKnownMarkers
#'
#' @param knownMarkers [knownMarkersDetected()] [tibble] grouped by cluster.
#' @param markdown Print Markdown headers.
#'
#' @return [writeLines()].
NULL



# Constructor ====
.plotKnownMarkers <- function(object, knownMarkers, markdown = TRUE) {
    if (nrow(knownMarkers) == 0L) {
        return(NULL)
    }
    cellTypes <- knownMarkers %>%
        pull("cellType") %>%
        unique
    pblapply(seq_along(cellTypes), function(a) {
        cellType <- cellTypes[[a]]
        if (isTRUE(markdown)) {
            mdHeader(cellType, level = 3L)
        }
        symbols <- knownMarkers %>%
            filter(.data[["cellType"]] == !!cellType) %>%
            pull("symbol") %>%
            unique %>%
            sort
        if (!is.null(symbols)) {
            plotClusters(object, symbols)
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
