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
#' @param object `data.frame` returned from [fetchTSNEExpressionData()].
#'
#' @return `ggplot`.
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
    color = scale_color_hue(),
    dark = TRUE,
    title = NULL
) {
    assert_is_data.frame(object)
    assert_is_character(axes)
    assert_is_subset(axes, colnames(object))
    assert_is_a_string(interestingGroups)
    assert_is_subset(interestingGroups, colnames(object))
    assert_is_a_bool(pointsAsNumbers)
    assert_is_a_number(pointSize)
    assert_is_a_number(pointAlpha)
    assert_is_a_bool(label)
    assert_is_a_number(labelSize)
    assertIsColorScaleDiscreteOrNULL(color)
    assert_is_a_bool(dark)
    assertIsAStringOrNULL(title)

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
            color = colorCol
        )
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
                    color = colorCol
                ),
                alpha = pointAlpha,
                size = pointSize
            )
    } else {
        p <- p +
            geom_point(
                alpha = pointAlpha,
                size = pointSize
            )
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
                    label = "ident"
                ),
                color = labelColor,
                size = labelSize,
                fontface = "bold"
            )
    }

    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    p
}
