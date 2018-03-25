#' Plot Top Markers
#'
#' @note The number of markers to plot is determined by the output of the
#' [topMarkers()] function. If you want to reduce the number of genes to plot,
#' simply reassign first using that function. If necessary, we can add support
#' for the number of genes to plot here in a future update.
#'
#' @name plotTopMarkers
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams plotMarkers
#' @inheritParams general
#' @param topMarkers `grouped_df` grouped by "`cluster`", returned by
#'   [topMarkers()].
#'
#' @return Show graphical output. Invisibly return `ggplot` plotlist.
#'
#' @examples
#' load(system.file(
#'     "extdata/all_markers_small.rda",
#'     package = "bcbioSingleCell"
#' ))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # seurat, grouped_df ====
#' plotTopMarkers(
#'     object = seurat_small,
#'     topMarkers = topMarkers(all_markers_small, n = 1L)
#' )
NULL



# Constructors =================================================================
.plotTopMarkers <- function(
    object,
    topMarkers,
    tsneColor = scale_color_viridis(discrete = FALSE),
    violinFill = scale_fill_viridis(discrete = TRUE),
    dotColor = scale_color_gradient(
        low = "lightgray",
        high = "purple"
    ),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    headerLevel = 2L
) {
    validObject(object)
    stopifnot(.isSanitizedMarkers(topMarkers))
    assertIsAHeaderLevel(headerLevel)

    clusters <- levels(topMarkers[["cluster"]])

    list <- pblapply(clusters, function(cluster) {
        genes <- topMarkers %>%
            .[.[["cluster"]] == cluster, , drop = FALSE] %>%
            pull("geneName")
        if (!length(genes)) {
            return(invisible())
        }
        if (length(genes) > 10L) {
            warn("Maximum of 10 genes per cluster is recommended")
        }
        markdownHeader(
            paste("Cluster", cluster),
            level = headerLevel,
            tabset = TRUE,
            asis = TRUE
        )
        subheaderLevel <- headerLevel + 1L
        plotMarkers(
            object = object,
            genes = genes,
            tsneColor = tsneColor,
            violinFill = violinFill,
            dotColor = dotColor,
            dark = dark,
            pointsAsNumbers = pointsAsNumbers,
            headerLevel = subheaderLevel
        )
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
        topMarkers = "grouped_df"
    ),
    .plotTopMarkers
)
