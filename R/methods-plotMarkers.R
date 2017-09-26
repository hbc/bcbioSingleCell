#' Plot Cell-Type Gene Markers
#'
#' Plot gene expression per cell in multiple formats:
#'
#' 1. t-SNE plot
#' 2. Violin plot
#' 3. Ridgeline plot
#' 4. Dot plot
#'
#' @rdname plotMarkers
#' @name plotMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#' @inheritParams plotMarkerTSNE
#' @param headerLevel Include a Markdown header for each gene.
#'
#' @return No value, only graphical output.
#'
#' @examples
#' \dontrun{
#' data(seurat, markers)
#' top <- topMarkers(markers)
#' genes <- top$symbol[1:4]
#' plotMarkers(seurat, genes = genes)
#' }
NULL



# Constructors ====
#' Plot Marker Seurat Constructor
#'
#' @param returnAsList Return the `gg` objects as a list instead of plotting
#'   with [cowplot::plot_grid()].
#'
#' @noRd
.plotMarkerSeurat <- function(object, gene, returnAsList = FALSE) {
    if (!is_string(gene)) {
        stop("gene must be a string", call. = FALSE)
    }

    lowExpressionCutoff <- 0.1

    # tSNE marker expression plot
    tsne <- plotMarkerTSNE(object,
                           genes = gene,
                           colorpoints = "expression")

    # Violin plot
    violin <- VlnPlot(
        object,
        features.plot = gene,
        do.return = TRUE,
        return.plotlist = TRUE) %>%
        .[[1]] %>%
        .[["data"]] %>%
        # Remove the low expression features
        dplyr::filter(.data[["feature"]] > lowExpressionCutoff) %>%
        ggplot(aes_(x = ~ident,
                    y = ~feature,
                    fill = ~ident)) +
        geom_violin(
            color = NA,
            scale = "width",
            adjust = 1,
            trim = TRUE) +
        scale_fill_viridis(discrete = TRUE) +
        theme(legend.position = "none")

    # Ridgeline (joy) plot
    ridges <- JoyPlot(
        object,
        features.plot = gene,
        cols.use = viridis(length(levels(object@ident))),
        do.return = TRUE,
        return.plotlist = TRUE) %>%
        .[[1]] %>%
        .[["data"]] %>%
        # Remove the low expression features
        dplyr::filter(.data[["feature"]] > lowExpressionCutoff) %>%
        ggplot(aes_(x = ~feature,
                    y = ~ident,
                    fill = ~ident)) +
        geom_density_ridges(color = NA, scale = 2) +
        scale_fill_viridis(discrete = TRUE) +
        theme(legend.position = "none")

    # Dot plot
    # `DotPlot()` is still plotting when `do.return = TRUE`
    dot <- DotPlot(
        object,
        genes.plot = gene,
        cols.use = c("white", "black"),
        do.return = TRUE,
        plot.legend = FALSE) %>%
        .[["data"]] %>%
        ggplot(aes_(x = ~genes.plot,
                     y = ~id)) +
        geom_point(aes_(size = ~pct.exp, color = ~avg.exp.scale)) +
        scale_radius(range = c(0, 6)) +
        scale_color_gradient(low = "white", high = "black") +
        labs(y = "ident") +
        theme(axis.title.x = element_blank(),
              legend.position = "none")

    if (isTRUE(returnAsList)) {
        list(tsne = tsne,
             violin = violin,
             ridges = ridges,
             dot = dot)
    } else {
        ggdraw() +
            # Coordinates are relative to lower left corner
            draw_plot(
                tsne,
                x = 0L, y = 0.3, width = 1, height = 0.7) +
            draw_plot(
                violin,
                x = 0, y = 0, width = 0.425, height = 0.3) +
            suppressMessages(draw_plot(
                ridges,
                x = 0.425, y = 0, width = 0.425, height = 0.3)) +
            draw_plot(
                dot,
                x = 0.85, y = 0, width = 0.15, height = 0.3)
    }
}



# Methods ====
#' @rdname plotMarkers
#' @export
setMethod("plotMarkers", "seurat", function(
    object,
    genes,
    headerLevel = 2L) {
    lapply(seq_along(genes), function(a) {
        gene <- genes[[a]]
        mdHeader(gene, level = headerLevel, asis = TRUE)
        .plotMarkerSeurat(object,
                          gene = gene,
                          returnAsList = FALSE) %>%
            show()
    }) %>%
        invisible()
})
