#' Plot t-SNE
#'
#' Generate a t-distributed stochastic neighbor embedding (t-SNE) plot.
#'
#' @rdname plotTSNE
#' @name plotTSNE
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @param interestingGroups Interesting group to use for plot colors.
#' @param pointsAsNumbers Plot the points as numbers (`TRUE`) or dots (`FALSE`).
#' @param pointSize Cell point size.
#' @param pointAlpha Alpha transparency level. Useful when there many cells in
#'   the dataset, and some cells can be masked.
#' @param label Overlay a cluster identitiy label on the plot.
#' @param labelSize Size of the text label.
#' @param color Color palette to use for points. Defaults to the viridis color
#'   palette. Use [ggplot2::scale_color_gradient()] to easily define your own
#'   low/high gradient.
#' @param dark Enable dark mode.
#' @param title Plot title.
#'
#' @return [ggplot].
#'
#' @examples
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # seurat
#' plotTSNE(seurat)
NULL



# Methods ======================================================================
#' @rdname plotTSNE
#' @importFrom ggplot2 scale_color_hue
#' @export
setMethod(
    "plotTSNE",
    signature("seurat"),
    function(
        object,
        interestingGroups = "ident",
        pointsAsNumbers = FALSE,
        pointSize = 0.5,
        pointAlpha = 0.8,
        label = TRUE,
        labelSize = 6L,
        color = ggplot2::scale_color_hue(),
        dark = TRUE,
        title = NULL) {
        tsne <- fetchTSNEData(object)
        .plotDR(
            tsne,
            axes = c(x = "tSNE1", y = "tSNE2"),
            interestingGroups = interestingGroups,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            color = color,
            dark = dark,
            title = title)
    })
