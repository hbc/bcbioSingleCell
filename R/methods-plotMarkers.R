#' Plot Gene Markers
#'
#' @rdname plotMarkers
#' @name plotMarkers
#'
#' @param symbols Character vector of gene marker symbols.
#' @param headerLevel Include a Markdown header for each gene.
#'
#' @return No value, only graphical output.
NULL



# Methods ====
#' @rdname plotMarkers
#' @export
setMethod("plotMarkers", "seurat", function(
    object,
    symbols,
    headerLevel = 2L) {
    lapply(seq_along(symbols), function(a) {
        message(symbols[[a]])
        mdHeader(symbols[[a]], level = headerLevel, asis = TRUE)
        # Joy plot
        JoyPlot(
            object,
            features.plot = symbols[[a]],
            cols.use = viridis(length(object@ident)),
            do.return = FALSE,
            x.lab.rot = TRUE) %>%
            show
        # Violin plot
        VlnPlot(
            object,
            features.plot = symbols[[a]],
            cols.use = viridis(length(object@ident)),
            do.return = FALSE,
            x.lab.rot = TRUE) %>%
            show
        # tSNE color plots
        # Dark theme enables greater contrast for marker visualization.
        # Otherwise use `rev(viridis(2))` for the colors. This will define
        # yellow as low and purple as high. The dark theme also shows
        FeaturePlot(
            object,
            features.plot = symbols[[a]],
            # Use viridis for better contrast than default colors
            # (1) low: purple; (2) high: yellow
            cols.use = viridis(2),
            do.return = FALSE,
            dark.theme = TRUE,
            no.legend = FALSE)
    }) %>%
        invisible
})
