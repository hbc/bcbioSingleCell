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



# Methods ======================================================================
setMethod(
    "plotFeature",
    signature("SingleCellExperiment"),
    function(
        object,
        features,
        reducedDim = c("TSNE", "UMAP"),
        color = getOption("bcbio.discrete.color", NULL),
        pointSize = getOption("bcbio.pointSize", 0.75),
        pointAlpha = getOption("bcbio.pointAlpha", 0.75),
        label = getOption("bcbio.label", TRUE),
        labelSize = getOption("bcbio.labelSize", 6L),
        dark = getOption("bcbio.dark", FALSE),
        grid = getOption("bcbio.grid", FALSE),
        legend = getOption("bcbio.legend", FALSE),
        aspectRatio = getOption("bcbio.aspectRatio", 1L)
    ) {
        assert_is_character(features)
        # Legacy support for `color = "auto"`
        if (identical(color, "auto")) {
            warning("Use `NULL` instead of `\"auto\"` for `color`")
            color <- NULL
        }
        assertIsColorScaleContinuousOrNULL(color)
        reducedDim <- match.arg(reducedDim)

        if (isTRUE(dark)) {
            fill <- "black"
        } else {
            fill <- "white"
        }

        data <- fetchReducedDimData(object, reducedDim = reducedDim)
        axes <- colnames(data)[seq_len(2L)]

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
                    x = !!sym("x"),
                    y = !!sym("y"),
                    color = !!sym(feature)
                )
            ) +
                geom_point(
                    alpha = pointAlpha,
                    size = pointSize
                ) +
                labs(
                    x = axes[[1L]],
                    y = axes[[2L]],
                    title = feature
                )

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
)
