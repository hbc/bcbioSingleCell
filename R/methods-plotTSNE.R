#' Plot t-SNE
#'
#' Generate a t-distributed stochastic neighbor embedding (t-SNE) plot.
#'
#' @rdname plotTSNE
#' @name plotTSNE
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param interestingGroups Interesting group to use for plot colors.
#' @param pointSize Cell point size.
#' @param label Overlay a cluster identitiy label on the plot.
#' @param labelSize Size of the text label.
#' @param dark Enable dark mode.
#'
#' @return [ggplot].
#'
#' @examples
#' \dontrun{
#' plotTSNE(seurat)
#' }
NULL



# Methods ====
#' @rdname plotTSNE
#' @importFrom ggplot2 scale_color_hue
#' @export
setMethod(
    "plotTSNE",
    signature("seurat"),
    function(
        object,
        interestingGroups = "ident",
        color = scale_color_hue(),
        pointSize = 1,
        label = TRUE,
        labelSize = 6,
        dark = TRUE) {
        tsne <- fetchTSNEData(object)
        .plotDimensionalityReduction(
            tsne,
            axes = c(x = "tSNE1", y = "tSNE2"),
            interestingGroups = interestingGroups,
            color = color,
            pointSize = pointSize,
            label = label,
            labelSize = labelSize,
            dark = dark)
    })
