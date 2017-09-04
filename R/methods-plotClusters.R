#' Plot Clusters
#'
#' @rdname plotClusters
#' @name plotClusters
#'
#' @param symbols Character vector of gene symbols.
#' @param headerLevel Include a Markdown header for each gene.
#'
#' @return No value, only graphical output.
NULL



# Methods ====
#' @rdname plotClusters
#' @export
setMethod("plotClusters", "seurat", function(
    object,
    symbols,
    headerLevel = 4L) {
    lapply(seq_along(symbols), function(a) {
        if (is.numeric(headerLevel)) {
            mdHeader(symbols[a], level = headerLevel, asis = TRUE)
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
            cols.use = c("grey", "purple"))
    }) %>%
        invisible
})
