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
    # TEMP Fix for Seurat issue #111
    if (any(grepl("\\-", symbols))) {
        warning("Seurat v2 fails on symbols starting with hyphen or number",
                call. = FALSE)
        symbols <- symbols %>%
            .[!grepl("\\-", .)]
        if (!length(symbols)) return(NULL)
    }
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
