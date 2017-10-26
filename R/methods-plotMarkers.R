#' Plot Cell-Type Gene Markers
#'
#' @description
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
#'
#' @param headerLevel Include a Markdown header for each gene.
#'
#' @return No value, only graphical output.
#'
#' @examples
#' \dontrun{
#' top <- topMarkers(markers)
#' genes <- top$symbol[1:4]
#' plotMarkers(seurat, genes = genes)
#' }
NULL



# Constructors ====
#' Plot Marker Seurat Constructor
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom cowplot draw_plot ggdraw
#' @importFrom dplyr filter
#' @importFrom ggplot2 aes_string element_blank geom_violin ggplot theme
#' @importFrom ggridges geom_density_ridges
#' @importFrom rlang is_string
#' @importFrom Seurat JoyPlot VlnPlot
#' @importFrom viridis scale_fill_viridis viridis
#'
#' @param returnAsList Return the `gg` objects as a list.
.plotMarkerSeurat <- function(
    object,
    gene,
    pointsAsNumbers = FALSE,
    returnAsList = FALSE) {
    if (!is_string(gene)) {
        stop("gene must be a string", call. = FALSE)
    }

    lowExpressionCutoff <- 0.1

    # tSNE marker expression plot
    tsne <- plotMarkerTSNE(
        object,
        genes = gene,
        colorPoints = "expression",
        pointsAsNumbers = pointsAsNumbers)

    # Violin plot
    violin <- VlnPlot(
        object,
        features.plot = gene,
        do.return = TRUE,
        return.plotlist = TRUE) %>%
        .[[1]] %>%
        .[["data"]] %>%
        # Remove the low expression features
        filter(.data[["feature"]] > lowExpressionCutoff) %>%
        ggplot(
            mapping = aes_string(
                x = "ident",
                y = "feature",
                fill = "ident")
        ) +
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
        filter(.data[["feature"]] > lowExpressionCutoff) %>%
        ggplot(
            mapping = aes_string(
                x = "feature",
                y = "ident",
                fill = "ident")
        ) +
        geom_density_ridges(color = NA, scale = 2) +
        scale_fill_viridis(discrete = TRUE) +
        theme(legend.position = "none")

    # Dot plot
    dot <- plotDot(object, genes = gene) +
        theme(axis.title.x = element_blank(),
              legend.position = "none")

    if (isTRUE(returnAsList)) {
        list(tsne = tsne,
             dot = dot,
             violin = violin,
             ridges = ridges)
    } else {
        ggdraw() +
            # Coordinates are relative to lower left corner
            draw_plot(
                tsne,
                x = 0, y = 0.25, width = 1, height = 0.75) +
            draw_plot(
                dot,
                x = 0, y = 0, width = 0.2, height = 0.25) +
            draw_plot(
                violin,
                x = 0.2, y = 0, width = 0.45, height = 0.25) +
            suppressMessages(draw_plot(
                ridges,
                x = 0.65, y = 0, width = 0.35, height = 0.25))
    }
}



# Methods ====
#' @rdname plotMarkers
#' @importFrom basejump mdHeader
#' @export
setMethod("plotMarkers", "seurat", function(
    object,
    genes,
    pointsAsNumbers = FALSE,
    headerLevel = 2) {
    lapply(seq_along(genes), function(a) {
        gene <- genes[[a]]
        mdHeader(gene, level = headerLevel, asis = TRUE)
        .plotMarkerSeurat(
            object,
            gene = gene,
            pointsAsNumbers = pointsAsNumbers,
            returnAsList = FALSE) %>%
            show()
    }) %>%
        invisible()
})
