#' Plot Known Markers Detected
#'
#' @rdname plotKnownMarkersDetected
#' @name plotKnownMarkersDetected
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inherit plotMarkers
#'
#' @param knownMarkersDetected [knownMarkersDetected()] return [tibble] grouped
#'   by cluster.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "knownMarkersDetected.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # seurat
#' # Let's plot the first 2 markers, as a quick example
#' plotKnownMarkersDetected(
#'     seurat,
#'     knownMarkersDetected = knownMarkersDetected[1:2, ])
NULL



# Constructors =================================================================
#' @importFrom basejump mdHeader
#' @importFrom dplyr pull
.plotKnownMarkersDetected <- function(
    object,
    knownMarkersDetected,
    color = scale_color_viridis(option = "inferno"),
    dark = TRUE,
    headerLevel = NULL) {
    if (!nrow(knownMarkersDetected)) return(NULL)
    cellTypes <- knownMarkersDetected %>%
        pull("cell") %>%
        unique()
    pblapply(seq_along(cellTypes), function(a) {
        cellType <- cellTypes[[a]]
        ensgene <- knownMarkersDetected %>%
            .[.[["cell"]] == cellType, ] %>%
            pull("ensgene") %>%
            unique()
        if (is.null(ensgene)) return(NULL)
        if (!is.null(headerLevel)) {
            mdHeader(
                cellType,
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE)
        }
        plotMarkers(
            object,
            genes = genes,
            format = "ensgene",
            color = color,
            dark = dark,
            headerLevel = headerLevel + 1)
        # `show()` is already declared in `plotMarkers()`
    }) %>%
        invisible()
}



# Methods ======================================================================
#' @rdname plotKnownMarkersDetected
#' @export
setMethod(
    "plotKnownMarkersDetected",
    signature(object = "seurat",
              knownMarkersDetected = "grouped_df"),
    .plotKnownMarkersDetected)
