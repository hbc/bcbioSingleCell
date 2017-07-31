#' Plot Top Markers
#'
#' @rdname plotTopMarkers
#' @name plotTopMarkers
#'
#' @param markers Top markers grouped [tibble] returned by [topMarkers()].
#' @param markdown Print Markdown headers.
#'
#' @return [ggplot].
#' @export
NULL



# Constructors ====
.plotTopMarkers <- function(object, markers, markdown = TRUE) {
    # Fix for gene symbol mismatch
    if ("gene" %in% colnames(markers)) {
        markers <- rename(markers, symbol = .data[["gene"]])
    }
    clusters <- markers[["cluster"]] %>% levels
    pblapply(seq_along(clusters), function(a) {
        cluster <- clusters[[a]]
        if (isTRUE(markdown)) {
            mdHeader(paste("Cluster", cluster), level = 3L)
        }
        symbols <- markers %>%
            # FIXME This isn't matching...
            .[.[["cluster"]] == cluster, ] %>%
            .[, "symbol"]
        if (is.null(symbols)) return(NULL)
        plotClusters(object, symbols)
    }) %>%
        invisible
}



# Methods ====
#' @rdname plotTopMarkers
#' @export
setMethod(
    "plotTopMarkers",
    signature(object = "seurat", markers = "grouped_df"),
    .plotTopMarkers)
