#' Plot Dimensionality Reduction
#'
#' Internal constructor supporting tSNE and PCA plots.
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom ggplot2 aes_string geom_point geom_text ggplot guide_legend
#'   guides labs scale_color_hue
#'
#' @param object [data.frame] returned from [fetchTSNEExpressionData()].
#'
#' @return [ggplot].
#' @noRd
.plotDimensionalityReduction <- function(
    object,
    axes,
    interestingGroups = "ident",
    color = scale_colour_hue(),
    pointSize = 1,
    labelSize = 6,
    dark = TRUE,
    label = TRUE) {
    if (interestingGroups == "ident") {
        # Seurat stores the ident from `FetchData()` as `object.ident`
        colorCol <- "ident"
    } else {
        colorCol <- interestingGroups
    }
    p <- ggplot(
        object,
        mapping = aes_string(
            x = axes[["x"]],
            y = axes[["y"]],
            color = colorCol)
    )
    # Put the dark theme call before the other ggplot aesthetics
    if (isTRUE(dark)) {
        p <- p + darkTheme()
    }
    p <- p +
        # Alpha transparency helps distinguish superimposed points
        geom_point(size = pointSize) +
        guides(color = guide_legend(title.position = "left", byrow = TRUE))
    if (!is.null(color)) {
        p <- p + color
    }
    if (interestingGroups == "ident") {
        # Change `ident` to `cluster` (more informative)
        p <- p + labs(color = "cluster")
    }
    if (isTRUE(label)) {
        if (isTRUE(dark)) {
            labelColor <- "white"
        } else {
            labelColor <- "black"
        }
        p <- p +
            geom_text(
                mapping = aes_string(
                    x = "centerX",
                    y = "centerY",
                    label = "ident"),
                color = labelColor,
                size = labelSize,
                # Hard-code this as bold (too hard to see otherwise)
                fontface = "bold")
    }
    p
}
