#' Plot Principal Component Analysis (PCA)
#'
#' @name plotPCA
#' @family Clustering Functions
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics plotPCA
#'
#' @inheritParams plotTSNE
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' # seurat ====
#' plotPCA(Seurat::pbmc_small)
NULL



# Methods ======================================================================
#' @rdname plotPCA
#' @export
setMethod(
    "plotPCA",
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
        title = NULL
    ) {
        pca <- fetchPCAData(object)
        .plotDR(
            pca,
            axes = c(x = "pc1", y = "pc2"),
            interestingGroups = interestingGroups,
            color = color,
            pointsAsNumbers = pointsAsNumbers,
            pointSize = pointSize,
            pointAlpha = pointAlpha,
            label = label,
            labelSize = labelSize,
            dark = dark,
            title = title
        )
    }
)
