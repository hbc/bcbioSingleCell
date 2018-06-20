# FIXME Rework this to use facet wrapping



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
            mapping = aes_string(
                x = dimCols[[1L]],
                y = dimCols[[2L]],
                color = feature
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

        if (isTRUE(dark)) {
            p <- p +
                theme_midnight(
                    aspect_ratio = aspectRatio,
                    grid = grid
                )
            if (is.null(color)) {
                color <- scale_colour_viridis(option = "plasma")
            }
        } else {
            p <- p +
                theme_paperwhite(
                    aspect_ratio = aspectRatio,
                    grid = grid
                )
            if (is.null(color)) {
                color <- scale_colour_viridis(begin = 1L, end = 0L)
            }
        }

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
