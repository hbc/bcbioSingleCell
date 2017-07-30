#' Plot Top Markers
#'
#' @rdname plotTopMarkers
#'
#' @param y Top markers grouped [tibble] returned by [topMarkers()].
#' @param markdown Print Markdown headers.
#'
#' @return [ggplot].
#' @export
setMethod(
    "plotTopMarkers",
    signature(x = "seurat", y = "grouped_df"),
    function(x, y, markdown = TRUE) {
        clusters <- y[["cluster"]] %>% levels
        pblapply(seq_along(clusters), function(a) {
            cluster <- clusters[[a]]
            if (isTRUE(markdown)) {
                mdHeader(paste("Cluster", cluster), level = 3L)
            }
            symbols <- y %>%
                .[.[["cluster"]] == cluster, ] %>%
                .[, "symbol"]
            if (!is.null(symbols)) {
                plotClusters(x, symbols)
            } else {
                NULL
            }
        }) %>%
            invisible
    })
