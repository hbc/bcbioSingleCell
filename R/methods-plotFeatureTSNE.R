# FIXME Fix return documentation.

#' Plot Feature t-SNE
#'
#' @name plotFeatureTSNE
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams plotTSNE
#' @param features Character vector of features (e.g. gene expression, PC
#'   scores, number of genes detected).
#' @param legend Show legends in paneled plots. Defaults to `FALSE` because
#'   typically these look too busy and the legends can get cut off.
#' @param return Return as "`plot`" or "`list`".
#'
#' @seealso [Seurat::FeaturePlot()].
#'
#' @return Depends on the return argument.
#'
#' @examples
#' # seurat ====
#' plotFeatureTSNE(pbmc_small, features = "PC1")
NULL



# Methods ======================================================================
#' @rdname plotFeatureTSNE
#' @importFrom Seurat FetchData
#' @importFrom basejump midnightTheme
#' @importFrom cowplot plot_grid
#' @importFrom ggplot2 aes_string geom_point ggplot scale_color_gradient theme
#' @export
setMethod(
    "plotFeatureTSNE",
    signature("seurat"),
    function(
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
        return = c("plot", "list")
    ) {
        return <- match.arg(return)

        tsne <- fetchTSNEData(object)
        data <- Seurat::FetchData(object, vars.all = features)
        assert_are_identical(rownames(tsne), rownames(data))
        assert_are_disjoint_sets(tsne, data)
        data <- cbind(tsne, data)
        # Remove duplicate columns
        data <- data[, unique(colnames(data))]
        plotlist <- lapply(features, function(feature) {
            p <- ggplot(
                data = data,
                mapping = aes_string(
                    x = "tSNE1",
                    y = "tSNE2",
                    color = feature
                )
            ) +
                geom_point(
                    alpha = pointAlpha,
                    size = pointSize
                ) +
                labs(title = feature)

            if (isTRUE(dark)) {
                p <- p + midnightTheme()
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

            if (!isTRUE(legend)) {
                p <- p + theme(legend.position = "none")
            }

            p
        })

        if (return == "plot") {
            plot_grid(plotlist = plotlist, labels = NULL)
        } else if (return == "list") {
            plotlist
        }
    }
)
