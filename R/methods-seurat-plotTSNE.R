#' Plot t-SNE
#'
#' Generate a t-distributed stochastic neighbor embedding (t-SNE) plot.
#'
#' @name plotTSNE
#' @family Clustering Functions
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' # seurat ====
#' plotTSNE(Seurat::pbmc_small)
NULL



# Methods ======================================================================
#' @rdname plotTSNE
#' @export
setMethod(
    "plotTSNE",
    signature("seurat"),
    function(
        object,
        interestingGroups = "ident",
        color = scale_color_hue(),
        pointsAsNumbers = FALSE,
        pointSize = 0.5,
        pointAlpha = 0.8,
        label = TRUE,
        labelSize = 6L,
        dark = TRUE,
        aspectRatio = 1L,
        title = NULL
    ) {
        tsne <- fetchTSNEData(object)
        .plotDR(
            tsne,
            axes = c(x = "tSNE1", y = "tSNE2"),
            interestingGroups = interestingGroups,
            color = color,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            dark = dark,
            aspectRatio = aspectRatio,
            title = title
        )
    }
)
