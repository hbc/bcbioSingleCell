#' Plot clusters
#'
#' @rdname plot_clusters
#' @author Michael Steinbaugh
#'
#' @param object Primary object.
#' @param symbols Character vector of gene symbols.
#'
#' @return No value, only graphical output.
#' @export
setMethod("plot_clusters", "seurat", function(object, symbols) {
    VlnPlot(object,
            features.plot = symbols,
            use.scaled = TRUE)
    FeaturePlot(object,
                features.plot = symbols,
                cols.use = c("grey", "blue"))
})



#' Plot known markers
#'
#' @rdname plot_known_markers
#'
#' @param x Primary object.
#' @param y [known_markers_detected()] [tibble] grouped by cluster.
#' @param markdown Print Markdown headers.
#'
#' @export
setMethod(
    "plot_known_markers",
    signature(x = "seurat", y = "grouped_df"),
    function(x, y, markdown = TRUE) {
        cell_type <- y %>%
            pull("cell_type") %>%
            unique
        pblapply(seq_along(cell_type), function(a) {
            if (isTRUE(markdown)) {
                # FIXME Add support for Markdown header level
                paste("###", cell_type[[a]]) %>%
                    print %>%
                    asis_output %>%
                    show
            }
            genes <- y %>%
                filter(.data[["cell_type"]] == cell_type[a]) %>%
                pull("symbol") %>%
                unique %>%
                sort
            plot_clusters(x, genes)
        }) %>%
            invisible
    })



#' Plot top markers
#'
#' @rdname plot_top_markers
#'
#' @param x Object.
#' @param y Top markers grouped [tibble] returned by [top_markers()].
#' @param markdown Print Markdown headers.
#'
#' @return [ggplot].
#' @export
setMethod(
    "plot_top_markers",
    signature(x = "seurat", y = "grouped_df"),
    function(x, y, markdown = TRUE) {
        cluster_levels <- y[["cluster"]] %>% levels
        pblapply(cluster_levels, function(a) {
            if (isTRUE(markdown)) {
                # FIXME Add support for Markdown header level
                paste("###", "Cluster", cluster_levels[a]) %>%
                    print %>%
                    asis_output %>%
                    show
            }
            genes <- y %>%
                filter(.data[["cluster"]] == a) %>%
                pull("symbol")
            plot_clusters(x, genes)
        }) %>%
            invisible
    })
