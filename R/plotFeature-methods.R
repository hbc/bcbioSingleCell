#' Plot Feature
#'
#' @name plotFeature
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param features `character`. Features to plot (e.g. gene expression, PC
#'   scores, number of genes detected).
#'
#' @seealso [Seurat::FeaturePlot()].
#'
#' @return `ggplot` or `list`.
#'
#' @examples
#' # SingleCellExperiment ====
#' object <- indrops_small
#' features <- c("nUMI", "nGene", "mitoRatio")
#'
#' # t-SNE
#' plotFeatureTSNE(object, features)
#'
#' # UMAP
#' plotFeatureUMAP(object, features)
NULL



# Constructors =================================================================
.plotFeatureReduction <- function(
    object,
    features,
    reduction = c("TSNE", "UMAP"),
    color = getOption("color", NULL),
    pointSize = getOption("pointSize", 0.75),
    pointAlpha = getOption("pointAlpha", 0.75),
    label = getOption("label", TRUE),
    labelSize = getOption("labelSize", 6L),
    dark = getOption("dark", FALSE),
    grid = getOption("grid", FALSE),
    legend = getOption("legend", FALSE),
    aspectRatio = getOption("aspectRatio", 1L)
) {
    assert_is_character(features)
    # Legacy support for `color = "auto"`
    if (identical(color, "auto")) {
        warning("Use `NULL` instead of `\"auto\"` for `color`")
        color <- NULL
    }
    assertIsColorScaleContinuousOrNULL(color)
    reduction <- match.arg(reduction)

    if (isTRUE(dark)) {
        fill <- "black"
    } else {
        fill <- "white"
    }

    # Get reduced dimension coordinates (t-SNE and UMAP are supported)
    if (reduction == "TSNE") {
        fxn <- fetchTSNEData
        dimCols <- c("tSNE_1", "tSNE_2")
    } else if (reduction == "UMAP") {
        fxn <- fetchUMAPData
        dimCols <- c("UMAP1", "UMAP2")
    }
    data <- fxn(object)

    # If the features are not defined, attempt to merge all reduced dims
    # information before stopping
    if (!all(features %in% colnames(data))) {
        reducedDimsData <- do.call(
            what = cbind,
            args = reducedDims(object)
        )
        stopifnot(identical(rownames(data), rownames(reducedDimsData)))
        data <- cbind(data, reducedDimsData) %>%
            # Ensure columns are unique
            .[, unique(colnames(.))]
    }
    assert_is_subset(features, colnames(data))

    plotlist <- lapply(features, function(feature) {
        p <- ggplot(
            data = data,
            mapping = aes(
                x = !!sym(dimCols[[1L]]),
                y = !!sym(dimCols[[2L]]),
                color = !!sym(feature)
            )
        ) +
            geom_point(
                alpha = pointAlpha,
                size = pointSize
            ) +
            labs(title = feature)

        if (isTRUE(label)) {
            if (isTRUE(dark)) {
                labelColor <- "white"
            } else {
                labelColor <- "black"
            }
            p <- p +
                geom_text(
                    mapping = aes(
                        x = !!sym("centerX"),
                        y = !!sym("centerY"),
                        label = !!sym("ident")
                    ),
                    color = labelColor,
                    size = labelSize,
                    fontface = "bold"
                )
        }

        if (isTRUE(dark)) {
            theme <- theme_midnight
            if (is.null(color)) {
                color <- darkMarkerColors
            }
        } else {
            theme <- theme_paperwhite
            if (is.null(color)) {
                color <- lightMarkerColors
            }
        }
        p <- p +
            theme(
                aspect_ratio = aspectRatio,
                grid = grid
            )

        if (is(color, "ScaleContinuous")) {
            p <- p + color
        }

        if (!isTRUE(legend)) {
            p <- p + guides(color = "none")
        }

        p
    })

    # Return ===================================================================
    if (length(features) > 1L) {
        plot_grid(plotlist = plotlist) +
            theme(
                plot.background = element_rect(
                    color = NA,
                    fill = fill
                )
            )
    } else {
        plotlist[[1L]]
    }
}



# Methods ======================================================================
#' @rdname plotFeature
#' @export
setMethod(
    "plotFeatureTSNE",
    signature("SingleCellExperiment"),
    function(
        object,
        features,
        color = getOption("color", NULL),
        pointSize = getOption("pointSize", 0.75),
        pointAlpha = getOption("pointAlpha", 0.75),
        label = getOption("label", TRUE),
        labelSize = getOption("labelSize", 6L),
        dark = getOption("dark", FALSE),
        grid = getOption("grid", FALSE),
        legend = getOption("legend", FALSE),
        aspectRatio = getOption("aspectRatio", 1L)
    ) {
        .plotFeatureReduction(
            object = object,
            features = features,
            reduction = "TSNE",
            color = color,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            dark = dark,
            grid = grid,
            legend = legend,
            aspectRatio = aspectRatio
        )
    }
)



#' @rdname plotFeature
#' @export
setMethod(
    "plotFeatureTSNE",
    signature("seurat"),
    getMethod("plotFeatureTSNE", "SingleCellExperiment")
)



#' @rdname plotFeature
#' @export
setMethod(
    "plotFeatureUMAP",
    signature("SingleCellExperiment"),
    function(
        object,
        features,
        color = getOption("color", NULL),
        pointSize = getOption("pointSize", 0.75),
        pointAlpha = getOption("pointAlpha", 0.75),
        label = getOption("label", TRUE),
        labelSize = getOption("labelSize", 6L),
        dark = getOption("dark", FALSE),
        grid = getOption("grid", FALSE),
        legend = getOption("legend", FALSE),
        aspectRatio = getOption("aspectRatio", 1L)
    ) {
        .plotFeatureReduction(
            object = object,
            features = features,
            reduction = "UMAP",
            color = color,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            dark = dark,
            grid = grid,
            legend = legend,
            aspectRatio = aspectRatio
        )
    }
)



#' @rdname plotFeature
#' @export
setMethod(
    "plotFeatureUMAP",
    signature("seurat"),
    getMethod("plotFeatureUMAP", "SingleCellExperiment")
)
