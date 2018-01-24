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
#' @param tsneColor Color palette to use for tSNE plot.
#' @param violinFill Color palette to use for violin plot.
#' @param dotColor Color palette to use for dot plot.
#' @param dark Plot the TSNE against a dark background.
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
#' @importFrom viridis scale_fill_viridis viridis
#'
#' @param returnAsList Return the `gg` objects as a list.
.plotMarker.seurat <- function(  # nolint
    object,
    gene,
    tsneColor = viridis::scale_color_viridis(),
    violinFill = viridis::scale_fill_viridis(discrete = TRUE),
    dotColor = ggplot2::scale_color_gradient(
        low = "lightgray",
        high = "purple"),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    title = NULL,
    return = "grid") {
    # Parameter integrity ======================================================
    if (!is_string(gene)) {
        abort("`gene` must be a string")
    }
    # return
    validReturn <- c("grid", "list")
    if (!return %in% validReturn) {
        abort(paste("`return` must contain:", toString(validReturn)))
    }

    # Plots ====================================================================
    tsne <- plotMarkerTSNE(
        object,
        genes = gene,
        format = "symbol",
        colorPoints = "expression",
        color = tsneColor,
        dark = dark,
        pointsAsNumbers = pointsAsNumbers,
        title = gene,
        subtitle = FALSE)
    violin <- plotViolin(
        object,
        genes = gene,
        format = "symbol",
        fill = violinFill)
    dot <- plotDot(
        object,
        genes = gene,
        format = "symbol",
        color = dotColor)

    # Return ===================================================================
    if (return == "grid") {
        dot <- dot +
            labs(x = "") +
            coord_flip() +
            theme(axis.text.y = element_text(angle = 90L, hjust = 0.5),
                  legend.position = "none")
        plot_grid(
            tsne,
            violin,
            dot,
            labels = NULL,
            ncol = 1L,
            nrow = 3L,
            rel_heights = c(1L, 0.3, 0.15)
        )
    } else if (return == "list") {
        list("tsne" = tsne,
             "dot" = dot,
             "violin" = violin)
    }
}



# Note the plural version here. This is the gene looping code.
.plotMarkers.seurat <- function(  # nolint
    object,
    genes,
    format = "symbol",
    tsneColor = viridis::scale_color_viridis(),
    violinFill = viridis::scale_fill_viridis(discrete = TRUE),
    dotColor = ggplot2::scale_color_gradient(
        low = "lightgray",
        high = "purple"),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    headerLevel = 2L,
    title = NULL) {
    .checkFormat(format)
    if (format == "ensgene") {
        genes <- .convertGenesToSymbols(object, genes = genes)
    }
    list <- lapply(seq_along(genes), function(a) {
        gene <- genes[[a]]
        # Skip and warn if gene is missing
        if (!gene %in% rownames(slot(object, "data"))) {
            return(warn(paste(gene, "missing")))
        }
        if (!is.null(headerLevel)) {
            mdHeader(gene, level = headerLevel, asis = TRUE)
        }
        p <- .plotMarker.seurat(
            object,
            gene = gene,
            tsneColor = tsneColor,
            violinFill = violinFill,
            dotColor = dotColor,
            dark = dark,
            pointsAsNumbers = pointsAsNumbers,
            title = title,
            return = "grid")
        show(p)
        p
    })
    invisible(list)
}



# Methods ======================================================================
#' @rdname plotMarkers
#' @importFrom bcbioBase mdHeader
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotMarkers",
    signature("seurat"),
    .plotMarkers.seurat)
