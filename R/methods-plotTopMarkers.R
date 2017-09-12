#' Plot Top Markers
#'
#' @rdname plotTopMarkers
#' @name plotTopMarkers
#' @inherit plotMarkers
#'
#' @param topMarkers Top markers grouped [tibble] returned by [topMarkers()].
NULL



# Constructors ====
.plotTopMarkers <- function(object, topMarkers, headerLevel = 2L) {
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
        if (is.null(symbols)) return(NULL)
        if (length(symbols) > 4L) {
            warning("Maximum of 4 genes per cluster is recommended")
            symbols <- symbols[[1L:4L]]
        }
        message(paste("Cluster", cluster))
        mdHeader(paste("Cluster", cluster),
                 level = headerLevel,
                 tabset = TRUE,
                 asis = TRUE)
        plotMarkers(object, symbols, headerLevel = headerLevel + 1L)
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
