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
#' @param returnAsList Return plots as a list.
#'
#' @seealso [Seurat::FeaturePlot()].
#'
#' @return No return, only graphical output.
#'
#' @examples
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # seurat
#' plotFeatureTSNE(seurat, features = "PC1")
NULL



# Methods ======================================================================
#' @rdname plotFeatureTSNE
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
        returnAsList = FALSE
    ) {
        data <- cbind(
            fetchTSNEData(object),
            FetchData(object, vars.all = features)
        )
        # Remove duplicate columns
        data <- data[, unique(colnames(data))]
        plotlist <- lapply(seq_along(features), function(a) {
            p <- ggplot(
                data,
                mapping = aes_string(
                    x = "tSNE1",
                    y = "tSNE2",
                    color = features[[a]]
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
                labs(title = features[[a]])
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
        if (isTRUE(returnAsList)) {
            plotlist
        } else {
            plot_grid(plotlist = plotlist, labels = NULL)
        }
    }
)
