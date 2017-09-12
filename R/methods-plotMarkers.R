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
            cols.use = viridis(length(symbols[[a]])),
            do.return = FALSE,
            features.plot = symbols[[a]],
            x.lab.rot = TRUE) %>%
            show
        # Violin plot
        VlnPlot(
            object,
            cols.use = viridis(length(symbols[[a]])),
            do.return = FALSE,
            features.plot = symbols[[a]],
            x.lab.rot = TRUE) %>%
            show
        # tSNE color plots
        FeaturePlot(
            object,
            do.return = FALSE,
            features.plot = symbols[[a]],
            # Use viridis for better contrast than default colors
            # (1) low: yellow; (2) high: purple
            cols.use = rev(viridis(2)))
    }) %>%
        invisible
})
