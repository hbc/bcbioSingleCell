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
#' @param label Label the clusters on the plot.
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
#' @export
setMethod(
    "plotTSNE",
    signature("seurat"),
    function(
        object,
        interestingGroups = "ident",
        label = TRUE,
        dark = TRUE) {
        tsne <- fetchTSNEData(object)
        .plotDimensionalityReduction(
            tsne,
            axes = c(x = "tSNE1", y = "tSNE2"),
            interestingGroups = interestingGroups,
            label = label,
            dark = dark)
    })
