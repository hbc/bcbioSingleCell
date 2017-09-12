#' Plot Gene Markers
#'
#' @rdname plotMarkers
#' @name plotMarkers
#'
#' @param symbols Character vector of gene marker symbols.
#' @param headerLevel Include a Markdown header for each gene.
#' @param combine Combine all markers into a single plot.
#'
#' @return No value, only graphical output.
NULL



# Constructors ====
.plotSeuratMarkers <- function(object, symbols, nCol = NULL) {
    # tSNE color plot
    # Dark theme enables greater contrast for marker visualization.
    # Otherwise use `rev(viridis(2))` for the colors. This will define
    # yellow as low and purple as high. The dark theme also shows
    FeaturePlot(
        object,
        features.plot = symbols,
        # Use viridis for better contrast than default colors
        # (1) low: purple; (2) high: yellow
        cols.use = viridis(2),
        do.return = FALSE,
        dark.theme = TRUE,
        nCol = nCol,
        no.legend = FALSE)

    # Violin plot
    VlnPlot(
        object,
        features.plot = symbols,
        cols.use = viridis(length(levels(object@ident))),
        do.return = FALSE,
        nCol = nCol,
        x.lab.rot = TRUE) %>%
        show

    # Joy plot
    JoyPlot(
        object,
        features.plot = symbols,
        cols.use = viridis(length(levels(object@ident))),
        do.return = FALSE,
        nCol = nCol) %>%
        show

    # Dot plot
    DotPlot(
        object,
        genes.plot = symbols,
        cols.use = viridis(2),
        plot.legend = TRUE,
        x.lab.rot = TRUE)
}



# Methods ====
#' @rdname plotMarkers
#' @export
setMethod("plotMarkers", "seurat", function(
    object,
    symbols,
    headerLevel = 2L,
    combine = TRUE) {
    if (isTRUE(combine)) {
        .plotSeuratMarkers(object, symbols, nCol = 2L)
    } else {
        lapply(seq_along(symbols), function(a) {
            mdHeader(symbols[[a]], level = headerLevel, asis = TRUE)
            .plotSeuratMarkers(object, symbols = symbols[[a]])
        }) %>%
            invisible
    }
})
