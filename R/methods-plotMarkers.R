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
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param tsneColor Color palette to use for tSNE plot.
#' @param violinFill Color palette to use for violin plot.
#' @param dotColor Color palette to use for dot plot.
#'
#' @return Show graphical output. Invisibly return `ggplot` plotlist.
#'
#' @examples
#' load(system.file(
#'     "extdata/all_markers_small.rda",
#'     package = "bcbioSingleCell"
#' ))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # seurat ===
#' top <- topMarkers(all_markers_small, n = 1L)
#' genes <- pull(top, "rowname")
#' plotMarkers(seurat_small, genes = genes)
NULL



# Constructors =================================================================
.plotMarker.seurat <- function(  # nolint
    object,
    gene,
    tsneColor = scale_color_viridis(discrete = FALSE),
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
    assertIsFillScaleDiscreteOrNULL(violinFill)
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



.plotMarkers.seurat <- function(  # nolint
    object,
    genes,
    tsneColor = scale_color_viridis(discrete = FALSE),
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
    # Passthrough: tsneColor, violinFill, dotColor, dark, pointsAsNumbers,
    # title
    assert_is_character(genes)
    assert_is_subset(genes, rownames(object))
    assertIsAHeaderLevel(headerLevel)

    list <- lapply(genes, function(gene) {
        markdownHeader(gene, level = headerLevel, asis = TRUE)
        p <- .plotMarker.seurat(
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
    .plotMarkers.seurat
)
