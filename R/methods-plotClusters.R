#' Plot Clusters
#'
#' @rdname plotClusters
#' @name plotClusters
#'
#' @param symbols Character vector of gene symbols.
#' @param markdown Include a Markdown header for each gene.
#'
#' @return No value, only graphical output.
NULL



# Methods ====
#' @rdname plotClusters
#' @export
setMethod("plotClusters", "seurat", function(
    object,
    symbols,
    markdown = TRUE) {
    sapply(seq_along(symbols), function(a) {
        if (isTRUE(markdown)) {
            mdHeader(symbols[a], level = 4L)
        }

        # Violin plots
        VlnPlot(
            object,
            do.return = FALSE,
            features.plot = symbols[a],
            x.lab.rot = TRUE) %>%
            show

        # tSNE color plots
        FeaturePlot(
            object,
            do.return = FALSE,
            features.plot = symbols[a],
            cols.use = c("grey", "blue"))
    }) %>%
        invisible
})
