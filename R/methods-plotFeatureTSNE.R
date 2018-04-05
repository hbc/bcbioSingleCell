#' Plot Feature t-SNE
#'
#' @name plotFeatureTSNE
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams plotTSNE
#' @inheritParams general
#' @param features Character vector of features (e.g. gene expression, PC
#'   scores, number of genes detected).
#' @param legend Show legends in paneled plots. Defaults to `FALSE` because
#'   typically these look too busy and the legends can get cut off.
#'
#' @seealso [Seurat::FeaturePlot()].
#'
#' @return `ggplot` or `list`.
#'
#' @examples
#' # seurat ====
#' plotFeatureTSNE(pbmc_small, features = "PC1")
NULL



# Methods ======================================================================
#' @rdname plotFeatureTSNE
#' @export
setMethod(
    "plotFeatureTSNE",
    signature("seurat"),
    function(  # nolint
        object,
        features,
        pointSize = 0.5,
        pointAlpha = 0.8,
        color = scale_color_gradient(
            low = "lightgray",
            high = "purple"
        ),
        dark = FALSE,
        label = TRUE,
        labelSize = 6L,
        legend = FALSE,
        return = c("ggplot", "list")
    ) {
        return <- match.arg(return)

        tsne <- fetchTSNEData(object)
        data <- FetchData(object, vars.all = features)
        assert_are_identical(rownames(tsne), rownames(data))
        assert_are_disjoint_sets(colnames(tsne), colnames(data))

        data <- cbind(tsne, data)

        plotlist <- lapply(features, function(feature) {
            p <- ggplot(
                data = data,
                mapping = aes_string(
                    x = "tSNE1",
                    y = "tSNE2",
                    color = feature
                )
            )

            if (isTRUE(dark)) {
                p <- p + midnightTheme()
            }

            p <- p +
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

            if (is(color, "ScaleContinuous")) {
                p <- p + color
            }

            if (!isTRUE(legend)) {
                p <- p + theme(legend.position = "none")
            }

            p
        })

        if (return == "ggplot") {
            plot_grid(plotlist = plotlist, labels = NULL)
        } else if (return == "list") {
            plotlist
        }
    }
)
