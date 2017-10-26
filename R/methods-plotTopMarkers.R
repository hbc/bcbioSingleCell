#' Plot Top Markers
#'
#' @note The number of markers to plot is determined by the output of the
#' [topMarkers()] function. If you want to reduce the number of genes to plot,
#' simply reassign first using that function. If necessary, we can add support
#' for the number of genes to plot here in a future update.
#'
#' @rdname plotTopMarkers
#' @name plotTopMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inherit plotMarkers
#'
#' @param topMarkers Top markers grouped [tibble] returned by [topMarkers()].
NULL



# Constructors ====
#' @importFrom dplyr pull rename
.plotTopMarkers <- function(
    object,
    topMarkers,
    pointsAsNumbers = FALSE,
    headerLevel = 2) {
    .checkSanitizedMarkers(topMarkers)
    clusters <- topMarkers[["cluster"]] %>%
        levels()
    pblapply(seq_along(clusters), function(a) {
        cluster <- clusters[[a]]
        genes <- topMarkers %>%
            .[.[["cluster"]] == cluster, ] %>%
            pull("symbol")
        if (is.null(genes)) {
            return(NULL)
        }
        if (length(genes) > 10) {
            warning("Maximum of 10 genes per cluster is recommended")
        }
        mdHeader(
            paste("Cluster", cluster),
            level = headerLevel,
            tabset = TRUE,
            asis = TRUE)
        plotMarkers(
            object,
            genes = genes,
            pointsAsNumbers = pointsAsNumbers,
            headerLevel = headerLevel + 1)
        # Don't show here, already defined in `plotMarkers()`
    }) %>%
        invisible()
}



# Methods ====
#' @rdname plotTopMarkers
#' @export
setMethod(
    "plotTopMarkers",
    signature(object = "seurat",
              topMarkers = "grouped_df"),
    .plotTopMarkers)
