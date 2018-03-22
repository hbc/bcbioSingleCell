#' Plot Marker t-SNE
#'
#' @rdname plotMarkerTSNE
#' @name plotMarkerTSNE
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams fetchTSNEExpressionData
#' @inheritParams plotTSNE
#'
#' @param expression Calculation to apply on the aggregate marker expression.
#'   Supports `mean` (default), `median`, or `sum`.
#' @param subtitle Include gene identifiers as subtitles.
#'
#' @return `ggplot`.
#'
#' @examples
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' genes <- counts(seurat) %>% rownames() %>% head()
#' print(genes)
#'
#' # seurat
#' plotMarkerTSNE(seurat, genes = genes)
#'
#' # data.frame
#' df <- fetchTSNEExpressionData(seurat, genes = genes)
#' plotMarkerTSNE(df, genes = genes)
NULL



# Constructors =================================================================
#' Plot Marker tSNE Constructor
#'
#' @keywords internal
#' @noRd
#'
#' @param object Marker gene expression `data.frame` returned from
#'   [fetchTSNEExpressionData()].
.plotMarkerTSNE <- function(
    object,
    genes,
    expression = "mean",
    pointsAsNumbers = FALSE,
    pointSize = 0.5,
    pointAlpha = 0.8,
    label = TRUE,
    labelSize = 6L,
    # FIXME Change color method to use function name?
    color = scale_color_viridis(),
    dark = TRUE,
    legend = TRUE,
    title = NULL,
    subtitle = TRUE
) {
    requiredCols <- c(
        "centerX",
        "centerY",
        "mean",
        "median",
        "ident",
        "sum",
        "tSNE1",
        "tSNE2")
    if (!all(requiredCols %in% colnames(object))) {
        abort(paste(
            "Required columns:", toString(requiredCols)
        ))
    }
    validExpression <- c("mean", "median", "sum")
    if (!expression %in% validExpression) {
        abort(paste(
            "`expression` must contain:",
            toString(validExpression)
        ))
    }

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
        object,
        mapping = aes_string(
            x = "tSNE1",
            y = "tSNE2",
            color = expression)
    )

    if (isTRUE(dark)) {
        p <- p + midnightTheme()
    }

    # Labels
    if (!is.character(title)) {
        title <- NULL
    }
    if (isTRUE(subtitle)) {
        subtitle <- genes
    } else {
        subtitle <- NULL
    }
    p <- p + labs(title = title, subtitle = subtitle)

    # Customize legend
    if (isTRUE(legend)) {
        p <- p +
            # Make the guide longer than normal, to improve appearance of values
            # containing a decimal point
            guides(
                color = guide_colorbar(
                    barwidth = 20L,
                    barheight = 1L,
                    direction = "horizontal")
            ) +
            theme(
                legend.justification = "center",
                legend.position = "bottom")
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
                    color = expression),
                alpha = pointAlpha,
                size = pointSize)
    } else {
        p <- p +
            geom_point(
                alpha = pointAlpha,
                size = pointSize)
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

    if (is(color, "ScaleContinuous")) {
        p <- p + color
    }

    p
}



# Methods ======================================================================
#' @rdname plotMarkerTSNE
#' @export
setMethod(
    "plotMarkerTSNE",
    signature("data.frame"),
    .plotMarkerTSNE)



#' @rdname plotMarkerTSNE
#' @export
setMethod(
    "plotMarkerTSNE",
    signature("seurat"),
    function(
        object,
        genes,
        expression = "mean",
        pointsAsNumbers = FALSE,
        pointSize = 0.5,
        pointAlpha = 0.8,
        label = TRUE,
        labelSize = 6L,
        color = scale_color_viridis(),
        dark = TRUE,
        title = NULL,
        subtitle = TRUE) {
        data <- fetchTSNEExpressionData(object, genes = genes)
        .plotMarkerTSNE(
            object = data,
            genes = genes,
            expression = expression,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            color = color,
            dark = dark,
            title = title,
            subtitle = subtitle
        )
    }
)
