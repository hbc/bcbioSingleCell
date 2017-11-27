#' Plot Principal Component Analysis (PCA)
#'
#' @rdname plotPCA
#' @name plotPCA
#' @family PCA Utilities
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics plotPCA
#'
#' @inherit plotTSNE
#'
#' @seealso [plotTSNE()].
NULL



# Methods ====
#' @rdname plotPCA
#' @importFrom ggplot2 scale_color_hue
#' @export
setMethod(
    "plotPCA",
    signature("seurat"),
    function(
        object,
        interestingGroups = "ident",
        pointsAsNumbers = FALSE,
        pointSize = 1,
        label = TRUE,
        labelSize = 6,
        color = scale_color_hue(),
        dark = TRUE,
        title = NULL) {
        pca <- fetchPCAData(object)
        .plotDimensionalityReduction(
            pca,
            axes = c(x = "pc1", y = "pc2"),
            interestingGroups = interestingGroups,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            label = label,
            labelSize = labelSize,
            color = color,
            dark = dark,
            title = title)
    })
