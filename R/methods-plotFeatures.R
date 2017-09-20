#' Plot Features of a Data Set
#'
#' @rdname plotFeatures
#' @name plotFeatures
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @param features Character vector of parameters supported by
#'   [Seurat::FetchData()] (e.g. nUMI, mitoRatio, PC).
#' @param headerLevel Include a Markdown header for each gene.
#' @param combine Combine all markers into a single plot.
#'
#' @return No value, only graphical output.
NULL



# Constructors ====
.plotFeaturesSeurat <- function(object, features, nCol = NULL) {
    # tSNE color plot
    # Dark theme enables greater contrast for marker visualization.
    # Otherwise use `rev(viridis(2))` for the colors. This will define
    # yellow as low and purple as high. The dark theme also shows
    tsne <- FeaturePlot(
        object,
        features.plot = features,
        # Use viridis for better contrast than default colors
        # (1) low: purple; (2) high: yellow
        cols.use = viridis(2),
        do.return = FALSE,
        dark.theme = TRUE,
        nCol = nCol,
        no.legend = FALSE)
    show(tsne)

    # Violin plot
    violin <- VlnPlot(
        object,
        features.plot = features,
        cols.use = viridis(length(levels(object@ident))),
        do.return = FALSE,
        nCol = nCol,
        x.lab.rot = TRUE)
    show(violin)

    # Joy plot
    joy <- JoyPlot(
        object,
        features.plot = features,
        cols.use = viridis(length(levels(object@ident))),
        do.return = FALSE,
        nCol = nCol)
    suppressMessages(show(joy))

    # Plots that are informative for only 2+ features
    if (length(features) > 1) {
        # Dot plot
        dot <- DotPlot(
            object,
            genes.plot = features,
            cols.use = viridis(2),
            do.return = FALSE,
            plot.legend = TRUE,
            x.lab.rot = TRUE)
        show(dot)
    }
}



# Methods ====
#' @rdname plotFeatures
#' @export
setMethod("plotFeatures", "seurat", function(
    object,
    features,
    headerLevel = 2L,
    combine = TRUE) {
    if (isTRUE(combine)) {
        .plotFeaturesSeurat(object, features, nCol = 2L)
    } else {
        lapply(seq_along(features), function(a) {
            mdHeader(features[[a]], level = headerLevel, asis = TRUE)
            .plotFeaturesSeurat(object, features[[a]])
        })
    }
    invisible()
})
