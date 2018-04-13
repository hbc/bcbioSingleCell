#' Plot Marker t-SNE
#'
#' @name plotMarkerTSNE
#' @family Clustering Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams fetchTSNEExpressionData
#' @inheritParams plotTSNE
#' @param expression Calculation to apply on the aggregate marker expression.
#' @param subtitle Include gene identifiers as subtitles.
#'
#' @return `ggplot`.
#'
#' @examples
#' # seurat ====
#' # Mitochondrial genes
#' mito <- grep("^MT\\.", rownames(counts(seurat_small)), value = TRUE)
#' head(mito)
#' plotMarkerTSNE(seurat_small, genes = mito)
NULL



# Methods ======================================================================
#' @rdname plotMarkerTSNE
#' @export
setMethod(
    "plotMarkerTSNE",
    signature("seurat"),
    .plotMarkerTSNE.seurat <- function(  # nolint
        object,
        genes,
        expression = c("mean", "median", "sum"),
        pointsAsNumbers = FALSE,
        pointSize = 0.5,
        pointAlpha = 0.8,
        label = TRUE,
        labelSize = 6L,
        color = scale_color_viridis(discrete = FALSE),
        dark = TRUE,
        legend = TRUE,
        title = NULL,
        subtitle = TRUE
    ) {
        assert_is_character(genes)
        expression <- match.arg(expression)
        assert_is_a_bool(pointsAsNumbers)
        assert_is_a_number(pointSize)
        assert_is_a_number(pointAlpha)
        assert_is_a_bool(label)
        assert_is_a_number(labelSize)
        assertIsColorScaleContinuousOrNULL(color)
        assert_is_a_bool(dark)
        assert_is_a_bool(legend)
        assertIsCharacterOrNULL(title)
        assert_is_a_bool(subtitle)

        data <- fetchTSNEExpressionData(object, genes = genes)
        requiredCols <- c(
            "centerX",
            "centerY",
            "mean",
            "median",
            "ident",
            "sum",
            "tSNE1",
            "tSNE2"
        )
        assert_is_subset(requiredCols, colnames(data))

        # Automatic subtitle containing list of marker genes
        if (isTRUE(subtitle)) {
            subtitle <- genes
            # Limit to the first 10 markers
            if (length(subtitle) > 10L) {
                subtitle <- c(subtitle[1L:10L], "...")
            }
            subtitle <- toString(subtitle)
        }

        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "tSNE1",
                y = "tSNE2",
                color = expression
            )
        )

        if (isTRUE(dark)) {
            p <- p + theme_midnight()
        } else {
            p <- p + theme_paperwhite()
        }

        # Labels
        if (isTRUE(subtitle)) {
            subtitle <- genes
        } else {
            subtitle <- NULL
        }
        p <- p + labs(title = title, subtitle = subtitle)

        # Customize legend
        if (isTRUE(legend)) {
            p <- p +
                # Make the guide longer than normal, to improve appearance of
                # values containing a decimal point
                guides(
                    color = guide_colorbar(
                        barwidth = 20L,
                        barheight = 1L,
                        direction = "horizontal"
                    )
                ) +
                theme(
                    legend.justification = "center",
                    legend.position = "bottom"
                )
        } else {
            p <- p + theme(legend.position = "none")
        }

        if (isTRUE(pointsAsNumbers)) {
            p <- p +
                geom_text(
                    mapping = aes_string(
                        x = "tSNE1",
                        y = "tSNE2",
                        label = "ident",
                        color = expression
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

        if (is(color, "ScaleContinuous")) {
            p <- p + color
        }

        p
    }
)
