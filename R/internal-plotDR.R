#' Plot Dimensionality Reduction
#'
#' Internal constructor supporting tSNE and PCA plots.
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @inheritParams plotTSNE
#'
#' @importFrom basejump midnightTheme
#' @importFrom ggplot2 aes_string geom_point geom_text ggplot guide_legend
#'   guides labs scale_color_hue
#'
#' @param object [data.frame] returned from [fetchTSNEExpressionData()].
#'
#' @return [ggplot].
#' @noRd
.plotDR <- function(
    object,
    axes,
    interestingGroups = "ident",
    pointsAsNumbers = FALSE,
    pointSize = 0.5,
    pointAlpha = 0.8,
    label = TRUE,
    labelSize = 6L,
    color = ggplot2::scale_color_hue(),
    dark = TRUE,
    title = NULL) {
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
        p <- p + midnightTheme()
    }
    if (isTRUE(pointsAsNumbers)) {
        p <- p +
            geom_text(
                mapping = aes_string(
                    x = axes[["x"]],
                    y = axes[["y"]],
                    label = "ident",
                    color = colorCol),
                alpha = pointAlpha,
                size = pointSize)
    } else {
        p <- p +
            geom_point(
                alpha = pointAlpha,
                size = pointSize)
    }
    if (interestingGroups == "ident") {
        # Present `ident` as `cluster` (more informative)
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
                fontface = "bold")
    }
    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }
    p
}
