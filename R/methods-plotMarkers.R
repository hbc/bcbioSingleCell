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
#' @param fill Color palette to use for violin plot fill. Defaults to viridis.
#' @param headerLevel Include a Markdown header for each gene.
#' @param title Plot title.
#'
#' @return No value, only graphical output.
#'
#' @examples
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "topMarkers.rda"),
#'     package = "bcbioSingleCell"))
#'
#' symbol <- topMarkers$symbol[[1]]
#' print(symbol)
#' ensgene <- topMarkers$ensgene[[1]]
#' print(ensgene)
#'
#' plotMarkers(seurat, genes = symbol, format = "symbol")
#' plotMarkers(seurat, genes = ensgene, format = "ensgene")
NULL



# Constructors =================================================================
#' Plot Marker Seurat Constructor
#'
#' @keywords internal
#' @noRd
#'
#' @importFrom cowplot plot_grid
#' @importFrom dplyr filter
#' @importFrom ggplot2 aes_string coord_flip element_blank geom_violin ggplot
#'   theme
#' @importFrom rlang is_string
#' @importFrom Seurat VlnPlot
#' @importFrom viridis scale_fill_viridis viridis
#'
#' @param returnAsList Return the `gg` objects as a list.
.plotMarkerSeurat <- function(
    object,
    gene,
    color = viridis::scale_color_viridis(),
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    title = NULL,
    returnAsList = FALSE) {
    if (!is_string(gene)) {
        stop("gene must be a string", call. = FALSE)
    }

    # tSNE marker expression plot
    tsne <- plotMarkerTSNE(
        object,
        genes = gene,
        format = "symbol",
        colorPoints = "expression",
        color = color,
        dark = dark,
        pointsAsNumbers = pointsAsNumbers,
        title = title)

    # Violin plot
    # FIXME Create our own `plotViolin()` function in a future update
    violin <- Seurat::VlnPlot(
        object,
        features.plot = gene,
        do.return = TRUE,
        return.plotlist = TRUE) %>%
        .[[1]] %>%
        .[["data"]] %>%
        # Remove the low expression features
        filter(.data[["feature"]] > 0.1) %>%
        ggplot(
            mapping = aes_string(
                x = "ident",
                y = "feature",
                fill = "ident")
        ) +
        geom_violin(
            color = "black",
            scale = "width",
            adjust = 1,
            trim = TRUE)
    if (is(fill, "ScaleDiscrete")) {
        violin <- violin + fill
    }


    # Dot plot
    # We're transposing the dot plot here to align vertically with the
    # violin and tSNE plots. The violin is preferable over the ridgeline
    # here because it works better horizontally.
    # FIXME Need to allow for user-defined color pass in here
    dot <- plotDot(object, genes = gene, format = "symbol")

    if (isTRUE(returnAsList)) {
        list(tsne = tsne,
             dot = dot,
             violin = violin)
    } else {
        # Customize the plots before preparing the grid
        violin <- violin +
            labs(y = "log expression") +
            theme(legend.position = "none")
        dot <- dot +
            labs(x = "") +
            coord_flip() +
            theme(axis.text.y = element_text(angle = 90, hjust = 0.5),
                  legend.position = "none")
        plot_grid(
            tsne,
            violin,
            dot,
            labels = NULL,
            ncol = 1,
            nrow = 3,
            rel_heights = c(1, 0.3, 0.15)
        )
    }
}



# Methods ======================================================================
#' @rdname plotMarkers
#' @importFrom basejump mdHeader
#' @importFrom viridis scale_color_viridis
#' @export
setMethod("plotMarkers", "seurat", function(
    object,
    genes,
    format = "symbol",
    color = viridis::scale_color_viridis(),
    fill = viridis::scale_fill_viridis(discrete = TRUE),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    headerLevel = NULL,
    title = NULL) {
    .checkFormat(format)
    if (format == "ensgene") {
        genes <- .convertGenesToSymbols(object, genes = genes)
    }
    return <- lapply(seq_along(genes), function(a) {
        gene <- genes[[a]]
        # Skip and warn if gene is missing
        if (!gene %in% rownames(slot(object, "data"))) {
            return(warning(paste(gene, "missing"), call. = FALSE))
        }
        if (!is.null(headerLevel)) {
            mdHeader(gene, level = headerLevel, asis = TRUE)
        }
        p <- .plotMarkerSeurat(
            object,
            gene = gene,
            color = color,
            fill = fill,
            dark = dark,
            pointsAsNumbers = pointsAsNumbers,
            title = title,
            returnAsList = FALSE)
        show(p)
        p
    })
    invisible(return)
})
