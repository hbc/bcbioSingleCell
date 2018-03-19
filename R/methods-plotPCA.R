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
#'
#' @examples
#' # seurat ====
#' plotPCA(pbmc_small)
NULL



# Methods ======================================================================
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
        pointSize = 0.5,
        pointAlpha = 0.8,
        label = TRUE,
        labelSize = 6L,
        color = scale_color_hue(),
        dark = TRUE,
        title = NULL
    ) {
        pca <- fetchPCAData(object)
        .plotDR(
            pca,
            axes = c(x = "pc1", y = "pc2"),
            interestingGroups = interestingGroups,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            color = color,
            dark = dark,
            title = title
        )
    }
)
