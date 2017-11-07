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
#' @export
setMethod(
    "plotPCA",
    signature("seurat"),
    function(
        object,
        interestingGroups = "ident",
        label = TRUE,
        dark = TRUE) {
        pca <- fetchPCAData(object)
        .plotDimensionalityReduction(
            pca,
            axes = c(x = "pc1", y = "pc2"),
            interestingGroups = interestingGroups,
            label = label,
            dark = dark)
    })
