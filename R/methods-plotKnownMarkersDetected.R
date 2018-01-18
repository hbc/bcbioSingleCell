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
#' @importFrom bcbioBase mdHeader
#' @importFrom dplyr pull
#' @importFrom stats na.omit
#' @importFrom viridis scale_color_viridis
.plotKnownMarkersDetected <- function(
    object,
    knownMarkersDetected,
    tsneColor = viridis::scale_color_viridis(),
    violinFill = viridis::scale_fill_viridis(discrete = TRUE),
    dotColor = ggplot2::scale_color_gradient(
        low = "lightgray",
        high = "purple"),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    headerLevel = 2) {
    if (!nrow(knownMarkersDetected)) return(NULL)
    cellTypes <- knownMarkersDetected %>%
        pull("cell") %>%
        unique() %>%
        na.omit()
    if (is.null(cellTypes)) return(NULL)
    list <- pblapply(seq_along(cellTypes), function(a) {
        cellType <- cellTypes[[a]]
        genes <- knownMarkersDetected %>%
            .[.[["cell"]] == cellType, ] %>%
            pull("symbol") %>%
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
            format = "symbol",
            tsneColor = tsneColor,
            violinFill = violinFill,
            dotColor = dotColor,
            dark = dark,
            pointsAsNumbers = pointsAsNumbers,
            headerLevel = subheaderLevel,
            title = cellType)
    })
    invisible(list)
}



# Methods ======================================================================
#' @rdname plotKnownMarkersDetected
#' @export
setMethod(
    "plotKnownMarkersDetected",
    signature(object = "seurat",
              knownMarkersDetected = "grouped_df"),
    .plotKnownMarkersDetected)
