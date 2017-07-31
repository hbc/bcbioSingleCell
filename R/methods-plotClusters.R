#' Plot Clusters
#'
#' @rdname plotClusters
#' @author Michael Steinbaugh
#'
#' @param symbols Character vector of gene symbols.
#'
#' @return No value, only graphical output.
#' @export
setMethod("plotClusters", "seurat", function(object, symbols) {
    VlnPlot(object,
            features.plot = symbols,
            use.scaled = TRUE)
    FeaturePlot(object,
                features.plot = symbols,
                cols.use = c("grey", "blue"))
})
