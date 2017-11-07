#' Plot Dimensionality Reduction
#'
#' Internal constructor supporting tSNE and PCA plots.
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @param object [data.frame] returned from [fetchTSNEExpressionData()].
#'
#' @return [ggplot].
#' @noRd
.plotDimensionalityReduction <- function(
    object,
    axes,
    interestingGroups = "ident",
    dark = TRUE,
    label = TRUE) {
    if (interestingGroups == "ident") {
        # Seurat stores the ident from `FetchData()` as `object.ident`
        color <- "ident"
    } else {
        color <- interestingGroups
    }
    p <- ggplot(
        object,
        mapping = aes_string(
            x = axes[["x"]],
            y = axes[["y"]],
            color = color)
    )
    # Put the dark theme call before the other ggplot aesthetics
    if (isTRUE(dark)) {
        p <- p + darkTheme()
    }
    p <- p +
        # Alpha transparency helps distinguish superimposed points
        geom_point() +
        guides(color = guide_legend(title.position = "left", byrow = TRUE))
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
                size = 6,
                fontface = "bold")
    }
    p
}
