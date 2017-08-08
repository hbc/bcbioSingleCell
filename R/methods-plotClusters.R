#' Plot Clusters
#'
#' @rdname plotClusters
#' @name plotClusters
#'
#' @param symbols Character vector of gene symbols.
#'
#' @return No value, only graphical output.
NULL



# Methods ====
#' @rdname plotClusters
#' @export
setMethod("plotClusters", "seurat", function(object, symbols) {
    VlnPlot(
        object,
        features.plot = symbols,
        nCol = 2L,
        x.lab.rot = TRUE)
    FeaturePlot(
        object,
        features.plot = symbols,
        cols.use = c("grey", "blue"),
        nCol = 2L)
})
