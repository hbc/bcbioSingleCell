#' Plot Top Markers
#'
#' @rdname plot_top_markers
#'
#' @param y Top markers grouped [tibble] returned by [top_markers()].
#' @param markdown Print Markdown headers.
#'
#' @return [ggplot].
#' @export
setMethod(
    "plot_top_markers",
    signature(x = "seurat", y = "grouped_df"),
    function(x, y, markdown = TRUE) {
        clusters <- y[["cluster"]] %>% levels
        pblapply(seq_along(clusters), function(a) {
            cluster <- clusters[[a]]
            if (isTRUE(markdown)) {
                mdHeader(paste("Cluster", cluster), level = 3L)
            }
            symbols <- y %>%
                filter(.data[["cluster"]] == !!cluster) %>%
                pull("symbol")
            if (!is.null(symbols)) {
                plot_clusters(x, symbols)
            } else {
                NULL
            }
        }) %>%
            invisible
    })
