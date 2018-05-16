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
#'
#' @return `ggplot`.
#'
#' @examples
#' # seurat ====
#' plotMarkerTSNE(seurat_small, genes = "COL1A1")
#' plotMarkerTSNE(seurat_small, genes = "COL1A1", dark = FALSE, grid = FALSE)
#'
#'
#' # Mitochondrial genes
#' mito <- grep("^MT\\.", rownames(counts(seurat_small)), value = TRUE)
#' print(sort(mito))
#' plotMarkerTSNE(seurat_small, genes = mito, title = "mitochondrial")
NULL



# Methods ======================================================================
#' @rdname plotMarkerTSNE
#' @export
setMethod(
    "plotMarkerTSNE",
    signature("seurat"),
    function(
        object,
        genes,
        expression = c("mean", "median", "sum"),
        color = "auto",
        pointsAsNumbers = FALSE,
        pointSize = 0.5,
        pointAlpha = 0.8,
        label = TRUE,
        labelSize = 6L,
        dark = FALSE,
        grid = TRUE,
        legend = FALSE,
        aspectRatio = 1L,
        title = TRUE
    ) {
        assert_is_character(genes)
        expression <- match.arg(expression)
        assert_is_a_bool(pointsAsNumbers)
        assert_is_a_number(pointSize)
        assert_is_a_number(pointAlpha)
        assert_is_a_bool(label)
        assert_is_a_number(labelSize)
        assert_is_a_bool(dark)
        assert_is_a_bool(legend)
        assert_is_any_of(title, c("character", "logical", "NULL"))
        if (is.character(title)) {
            assert_is_a_string(title)
        }

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

        p <- ggplot(
            data = data,
            mapping = aes_string(
                x = "tSNE1",
                y = "tSNE2",
                color = expression
            )
        )

        # Titles
        subtitle <- NULL
        if (isTRUE(title)) {
            if (is_a_string(genes)) {
                title <- genes
            } else {
                title <- NULL
                subtitle <- genes
                # Limit to the first 5 markers
                if (length(subtitle) > 5L) {
                    subtitle <- c(subtitle[1L:5L], "...")
                }
                subtitle <- toString(subtitle)
            }
        } else if (identical(title, FALSE)) {
            title <- NULL
        }
        p <- p + labs(title = title, subtitle = subtitle)

        # Customize legend
        if (isTRUE(legend)) {
            # Make the guide longer than normal, to improve appearance of values
            # containing a decimal point
            p <- p +
                guides(
                    color = guide_colorbar(
                        barwidth = 20L,
                        barheight = 1L,
                        direction = "horizontal"
                    )
                )
        } else {
            p <- p + guides(color = "none")
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

        # Color palette
        if (isTRUE(dark)) {
            p <- p +
                theme_midnight(
                    aspect_ratio = aspectRatio,
                    grid = grid
                )
            if (color == "auto") {
                color <- scale_color_viridis(
                    option = "plasma",
                    discrete = FALSE
                )
            }
        } else {
            p <- p +
                theme_paperwhite(
                    aspect_ratio = aspectRatio,
                    grid = grid
                )
            if (color == "auto") {
                color <- scale_color_gradient(
                    low = "gray90",
                    high = "black"
                )
            }
        }

        if (is(color, "ScaleContinuous")) {
            p <- p + color
        }

        p
    }
)
