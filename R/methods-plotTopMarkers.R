#' Plot Top Markers
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
    # Fix for gene symbol mismatch
    if ("gene" %in% colnames(topMarkers)) {
        topMarkers <- rename(topMarkers, symbol = .data[["gene"]])
    }
    clusters <- topMarkers[["cluster"]] %>% levels
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
        mdHeader(paste("Cluster", cluster),
                 level = headerLevel,
                 tabset = TRUE,
                 asis = TRUE)
        plotMarkers(object,
                    genes = genes,
                    pointsAsNumbers = pointsAsNumbers,
                    headerLevel = headerLevel + 1)
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
