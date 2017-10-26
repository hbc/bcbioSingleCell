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
#' @importFrom basejump mdHeader
#' @importFrom dplyr filter pull
.plotKnownMarkers <- function(
    object,
    knownMarkers,
    color = scale_color_viridis(option = "inferno"),
    dark = TRUE,
    headerLevel = 2) {
    if (nrow(knownMarkers) == 0) {
        return(NULL)
    }
    cellTypes <- knownMarkers %>%
        pull("cell") %>%
        unique()
    pblapply(seq_along(cellTypes), function(a) {
        cellType <- cellTypes[[a]]
        genes <- knownMarkers %>%
            filter(.data[["cell"]] == !!cellType) %>%
            pull("symbol") %>%
            unique() %>%
            sort()
        if (!is.null(genes)) {
            mdHeader(cellType, level = headerLevel, tabset = TRUE, asis = TRUE)
            plotMarkers(
                object,
                genes = genes,
                color = color,
                dark = dark,
                headerLevel = headerLevel + 1)
            # Show is already declared in `plotMarkers()`
        }
    }) %>%
        invisible()
}



# Methods ====
#' @rdname plotKnownMarkers
#' @export
setMethod(
    "plotKnownMarkers",
    signature(object = "seurat", knownMarkers = "grouped_df"),
    .plotKnownMarkers)
