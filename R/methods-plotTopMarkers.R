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
.plotTopMarkers <- function(
    object,
    topMarkers,
    headerLevel = 2,
    combine = FALSE) {
    # Fix for gene symbol mismatch
    if ("gene" %in% colnames(topMarkers)) {
        topMarkers <- rename(topMarkers, symbol = .data[["gene"]])
    }
    clusters <- topMarkers[["cluster"]] %>% levels
    pblapply(seq_along(clusters), function(a) {
        cluster <- clusters[[a]]
        symbols <- topMarkers %>%
            .[.[["cluster"]] == cluster, ] %>%
            pull("symbol")
        if (is.null(symbols)) {
            return(NULL)
        }
        if (length(symbols) > 10) {
            warning("Maximum of 10 genes per cluster is recommended")
        }
        mdHeader(paste("Cluster", cluster),
                 level = headerLevel,
                 tabset = TRUE,
                 asis = TRUE)
        plotMarkers(object,
                    symbols = symbols,
                    headerLevel = headerLevel + 1,
                    combine = combine)
    }) %>%
        invisible
}



# Methods ====
#' @rdname plotTopMarkers
#' @export
setMethod(
    "plotTopMarkers",
    signature(object = "seurat", topMarkers = "grouped_df"),
    .plotTopMarkers)
