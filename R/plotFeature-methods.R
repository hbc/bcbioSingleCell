#' Plot Feature
#'
#' @name plotFeature
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param features Character vector of features (e.g. gene expression, PC
#'   scores, number of genes detected).
#'
#' @seealso [Seurat::FeaturePlot()].
#'
#' @return `ggplot` or `list`.
#'
#' @examples
#' # seurat ====
#' object <- seurat_small
#' features <- c("nUMI", "nGene", "PC1", "PC2")
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
    color = NULL,
    pointSize = getOption("pointSize", 0.75),
    pointAlpha = getOption("pointAlpha", 0.75),
    label = TRUE,
    labelSize = getOption("labelSize", 6L),
    dark = FALSE,
    grid = FALSE,
    legend = FALSE,
    aspectRatio = 1L
) {
    assert_is_character(features)
    # Legacy support for `color = "auto"`
    if (identical(color, "auto")) {
        color <- NULL
    }
    assertIsColorScaleContinuousOrNULL(color)
    reduction <- match.arg(reduction)

    if (isTRUE(dark)) {
        fill <- "black"
    } else {
        fill <- "white"
    }

    if (reduction == "TSNE") {
        fxn <- fetchTSNEData
        dimCols <- c("tSNE_1", "tSNE_2")
    } else if (reduction == "UMAP") {
        fxn <- fetchUMAPData
        dimCols <- c("UMAP1", "UMAP2")
    }
    reducedDims <- fxn(object)

    featureData <- FetchData(object, vars.all = features)

    # Columns from `FetchData` take priority, if there is overlap
    if (length(intersect(colnames(reducedDims), colnames(featureData)))) {
        reducedDims <- reducedDims %>%
            .[, setdiff(colnames(.), colnames(featureData))]
    }
    assert_are_identical(rownames(reducedDims), rownames(featureData))
    assert_are_disjoint_sets(colnames(reducedDims), colnames(featureData))
    data <- cbind(reducedDims, featureData)

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
                color <- scale_colour_viridis(option = "plasma")
            }
        } else {
            theme <- theme_paperwhite
            # FIXME Just use gray/red here by default?
            if (is.null(color)) {
                color <- scale_colour_viridis(begin = 1L, end = 0L)
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
        color = NULL,
        pointSize = getOption("pointSize", 0.75),
        pointAlpha = getOption("pointAlpha", 0.75),
        label = TRUE,
        labelSize = getOption("labelSize", 6L),
        dark = FALSE,
        grid = FALSE,
        legend = FALSE,
        aspectRatio = 1L
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
        color = NULL,
        pointSize = getOption("pointSize", 0.75),
        pointAlpha = getOption("pointAlpha", 0.75),
        label = TRUE,
        labelSize = getOption("labelSize", 6L),
        dark = FALSE,
        grid = FALSE,
        legend = FALSE,
        aspectRatio = 1L
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
