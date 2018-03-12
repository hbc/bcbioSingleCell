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
#' @param topMarkers Top markers [tibble] grouped by cluster, returned by
#'   [topMarkers()].
#'
#' @examples
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/topMarkers.rda", package = "bcbioSingleCell"))
#'
#' # seurat, grouped_df
#' # Let's plot the top 2 markers from cluster 0, as a quick example
#' plotTopMarkers(seurat, topMarkers = topMarkers[1:2, ])
NULL



# Constructors =================================================================
#' @importFrom basejump markdownHeader
#' @importFrom dplyr rename
#' @importFrom pbapply pblapply
.plotTopMarkers <- function(
    object,
    topMarkers,
    tsneColor = viridis::scale_color_viridis(),
    violinFill = viridis::scale_fill_viridis(discrete = TRUE),
    dotColor = ggplot2::scale_color_gradient(
        low = "lightgray",
        high = "purple"),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    headerLevel = 2L) {
    .checkSanitizedMarkers(topMarkers)
    clusters <- levels(topMarkers[["cluster"]])
    list <- pblapply(seq_along(clusters), function(a) {
        cluster <- clusters[[a]]
        # We're matching against the `geneName` column here
        genes <- topMarkers %>%
            as.data.frame() %>%
            .[.[["cluster"]] == cluster, "geneName", drop = TRUE]
        if (!length(genes)) return(invisible())
        if (length(genes) > 10L) {
            warn("Maximum of 10 genes per cluster is recommended")
        }
        if (!is.null(headerLevel)) {
            markdownHeader(
                paste("Cluster", cluster),
                level = headerLevel,
                tabset = TRUE,
                asis = TRUE)
            subheaderLevel <- headerLevel + 1L
        } else {
            subheaderLevel <- NULL
        }
        plotMarkers(
            object,
            genes = genes,
            tsneColor = tsneColor,
            violinFill = violinFill,
            dotColor = dotColor,
            dark = dark,
            pointsAsNumbers = pointsAsNumbers,
            headerLevel = subheaderLevel)
    })
    invisible(list)
}



# Methods ======================================================================
#' @rdname plotTopMarkers
#' @export
setMethod(
    "plotTopMarkers",
    signature(
        object = "seurat",
        topMarkers = "grouped_df"),
    .plotTopMarkers)
