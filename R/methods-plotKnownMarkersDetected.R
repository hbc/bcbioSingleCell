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
#' @importFrom stats na.omit
#' @importFrom viridis scale_color_viridis
.plotKnownMarkersDetected <- function(
    object,
    knownMarkersDetected,
    color = viridis::scale_color_viridis(),
    dark = TRUE,
    headerLevel = NULL) {
    if (!nrow(knownMarkersDetected)) return(NULL)
    cellTypes <- knownMarkersDetected %>%
        pull("cell") %>%
        unique() %>%
        na.omit()
    if (is.null(cellTypes)) return(NULL)
    return <- pblapply(seq_along(cellTypes), function(a) {
        cellType <- cellTypes[[a]]
        genes <- knownMarkersDetected %>%
            .[.[["cell"]] == cellType, ] %>%
            pull("ensgene") %>%
            unique() %>%
            na.omit()
        if (is.null(genes)) return(NULL)
        if (!is.null(headerLevel)) {
            mdHeader(
                cellType,
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE)
            subheaderLevel <- headerLevel + 1
        } else {
            subheaderLevel <- NULL
        }
        plotMarkers(
            object,
            genes = genes,
            format = "ensgene",
            color = color,
            dark = dark,
            headerLevel = subheaderLevel,
            title = cellType)
    })
    invisible(return)
}



# Methods ======================================================================
#' @rdname plotKnownMarkersDetected
#' @export
setMethod(
    "plotKnownMarkersDetected",
    signature(object = "seurat",
              knownMarkersDetected = "grouped_df"),
    .plotKnownMarkersDetected)
