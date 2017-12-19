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
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "topMarkers.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # seurat, grouped_df
#' # Let's plot the top 2 markers from cluster 0, as a quick example
#' plotTopMarkers(seurat, topMarkers[1:2, ])
NULL



# Constructors =================================================================
#' @importFrom dplyr rename
#' @importFrom pbapply pblapply
.plotTopMarkers <- function(
    object,
    topMarkers,
    pointsAsNumbers = FALSE,
    headerLevel = 2) {
    .checkSanitizedMarkers(topMarkers)
    clusters <- levels(topMarkers[["cluster"]])
    return <- pblapply(seq_along(clusters), function(a) {
        cluster <- clusters[[a]]
        # We're matching against the `symbol` column here
        genes <- topMarkers %>%
            as.data.frame() %>%
            .[.[["cluster"]] == cluster, "symbol", drop = TRUE]
        if (is.null(genes)) return(NULL)
        if (length(genes) > 10) {
            warning(
                "Maximum of 10 genes per cluster is recommended",
                call. = FALSE)
        }
        if (!is.null(headerLevel)) {
            mdHeader(
                paste("Cluster", cluster),
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
            pointsAsNumbers = pointsAsNumbers,
            headerLevel = subheaderLevel)
    })
    invisible(return)
}



# Methods ======================================================================
#' @rdname plotTopMarkers
#' @export
setMethod(
    "plotTopMarkers",
    signature(object = "seurat",
              topMarkers = "grouped_df"),
    .plotTopMarkers)
