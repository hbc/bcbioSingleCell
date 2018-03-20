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
#' @name plotMarkers
#' @family Clustering Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams plotMarkerTSNE
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
#' # seurat ===
#' top <- topMarkers(all_markers_small, n = 1L)
#' genes <- pull(top, "geneName")
#' plotMarkers(seurat_small, genes = genes)
NULL



# Constructors =================================================================
# FIXME Export `plotMarker()`?
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 aes_string coord_flip element_blank geom_violin ggplot
#'   theme
.plotMarker <- function(  # nolint
    object,
    gene,
    tsneColor = scale_color_viridis(),
    violinFill = scale_fill_viridis(discrete = TRUE),
    dotColor = scale_color_gradient(
        low = "lightgray",
        high = "purple"
    ),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    title = NULL,
    return = c("grid", "list")
) {
    assert_is_a_string(gene)
    assertIsColorScaleContinuousOrNULL(tsneColor)
    assertIsColorScaleDiscreteOrNULL(violinFill)
    assertIsColorScaleContinuousOrNULL(dotColor)
    assert_is_a_bool(dark)
    assert_is_a_bool(pointsAsNumbers)
    assertIsAStringOrNULL(title)
    return <- match.arg(return)

    # Plots ====================================================================
    tsne <- plotMarkerTSNE(
        object,
        genes = gene,
        expression = "sum",
        color = tsneColor,
        dark = dark,
        pointsAsNumbers = pointsAsNumbers,
        title = gene,
        subtitle = FALSE
    )
    violin <- plotViolin(
        object,
        genes = gene,
        fill = violinFill
    )
    dot <- plotDot(
        object,
        genes = gene,
        color = dotColor
    )

    # Return ===================================================================
    if (return == "grid") {
        dot <- dot +
            labs(x = "") +
            coord_flip() +
            theme(
                axis.text.y = element_text(angle = 90L, hjust = 0.5),
                legend.position = "none"
            )
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
        list(
            "tsne" = tsne,
            "dot" = dot,
            "violin" = violin
        )
    }
}



#' @importFrom basejump markdownHeader
.plotMarkers <- function(  # nolint
    object,
    genes,
    tsneColor = scale_color_viridis(),
    violinFill = scale_fill_viridis(discrete = TRUE),
    dotColor = scale_color_gradient(
        low = "lightgray",
        high = "purple"
    ),
    dark = TRUE,
    pointsAsNumbers = FALSE,
    headerLevel = 2L,
    title = NULL
) {
    assertIsAHeaderLevel(headerLevel)

    # Abort on missing genes
    # FIXME Add method support for `rownames()` on seurat object
    assert_is_subset(genes, rownames(counts(object)))

    list <- lapply(genes, function(gene) {
        markdownHeader(gene, level = headerLevel, asis = TRUE)
        p <- .plotMarker(
            object = object,
            gene = gene,
            tsneColor = tsneColor,
            violinFill = violinFill,
            dotColor = dotColor,
            dark = dark,
            pointsAsNumbers = pointsAsNumbers,
            title = title,
            return = "grid"
        )
        show(p)
        invisible(p)
    })

    invisible(list)
}



# Methods ======================================================================
#' @rdname plotMarkers
#' @export
setMethod(
    "plotMarkers",
    signature("seurat"),
    .plotMarkers
)
