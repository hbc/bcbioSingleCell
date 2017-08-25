#' Plot Top Markers
#'
#' @rdname plotTopMarkers
#' @name plotTopMarkers
#'
#' @param topMarkers Top markers grouped [tibble] returned by [topMarkers()].
#' @param markdown Print Markdown headers.
#'
#' @return [ggplot].
NULL



# Constructors ====
.plotTopMarkers <- function(object, topMarkers, markdown = TRUE) {
    # Fix for gene symbol mismatch
    if ("gene" %in% colnames(topMarkers)) {
        topMarkers <- rename(topMarkers, symbol = .data[["gene"]])
    }
    clusters <- topMarkers[["cluster"]] %>% levels
    pblapply(seq_along(clusters), function(a) {
        cluster <- clusters[[a]]
        if (isTRUE(markdown)) {
            mdHeader(paste("Cluster", cluster), level = 3L, tabset = TRUE)
        }
        symbols <- topMarkers %>%
            .[.[["cluster"]] == cluster, ] %>%
            pull("symbol")
        if (is.null(symbols)) return(NULL)
        if (length(symbols) > 4L) {
            warning("Maximum of 4 genes per cluster is recommended")
            symbols <- symbols[[1L:4L]]
        }
        plotClusters(object, symbols)
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
