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
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # seurat
#' plotPCA(seurat)
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
        pointSize = 1L,
        label = TRUE,
        labelSize = 6L,
        color = ggplot2::scale_color_hue(),
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
