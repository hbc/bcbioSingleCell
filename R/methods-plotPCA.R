#' Plot Principal Component Analysis (PCA)
#'
#' @rdname plotPCA
#' @name plotPCA
#' @family PCA Utilities
#' @author Michael Steinbaugh
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
        interestingGroup = "ident",
        label = TRUE) {
        pca <- fetchPCAData(object)
        .plotDim(
            pca,
            axes = c(x = "pc1", y = "pc2"),
            interestingGroup = interestingGroup,
            label = label)
    })
