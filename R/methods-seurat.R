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
        cell_types <- y %>%
            pull("cell_type") %>%
            unique
        pblapply(seq_along(cell_types), function(a) {
            cell_type <- cell_types[[a]]
            if (isTRUE(markdown)) {
                # FIXME Add support for Markdown header level
                paste("###", cell_type) %>%
                    asis_output %>%
                    show
            }
            symbols <- y %>%
                filter(.data[["cell_type"]] == !!cell_type) %>%
                pull("symbol") %>%
                unique %>%
                sort
            if (is.null(symbols)) {
                plot_clusters(x, symbols)
            } else {
                NULL
            }
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
        clusters <- y[["cluster"]] %>% levels
        pblapply(seq_along(clusters), function(a) {
            cluster <- clusters[[a]]
            if (isTRUE(markdown)) {
                paste("###", "Cluster", cluster) %>%
                    asis_output %>%
                    show
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
