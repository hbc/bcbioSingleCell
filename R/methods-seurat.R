#' Plot clusters
#'
#' @rdname plot_clusters
#' @author Michael Steinbaugh
#'
#' @param object Primary object.
#' @param genes Character vector of gene symbols.
#'
#' @return No value, only graphical output.
#' @export
setMethod("plot_clusters", "seurat", function(object, genes) {
    VlnPlot(
        seurat,
        features.plot = genes,
        use.scaled = TRUE)
    FeaturePlot(
        seurat,
        features.plot = genes,
        cols.use = c("grey", "blue"))
})



#' Plot known markers
#'
#' @rdname plot_known_markers
#'
#' @param object Primary object.
#' @param known_markers_detected Known markers detected [data.frame].
#'
#' @export
setMethod("plot_known_markers", "seurat", function(
    object, known_markers_detected) {
    cell_type <- known_markers_detected %>% pull("cell_type") %>% unique
    pblapply(seq_along(cell_type), function(a) {
        genes <- known_markers_detected %>%
            filter(.data[["cell_type"]] == cell_type[a]) %>%
            pull("gene") %>%
            unique %>%
            sort
        plot_clusters(seurat, genes)
    }) %>%
        invisible
})



#' Plot top markers
#'
#' @rdname plot_top_markers
#'
#' @param object Primary object.
#' @param top_markers [data.frame] returned by [top_markers()].
#'
#' @return [ggplot].
#' @export
setMethod("plot_top_markers", "seurat", function(
    object, top_markers) {
    pblapply(levels(top_markers[["cluster"]]), function(a) {
        genes <- top_markers %>%
            filter(.data[["cluster"]] == a) %>%
            pull("gene")
        plot_clusters(seurat, genes)
    }) %>%
        invisible
})
