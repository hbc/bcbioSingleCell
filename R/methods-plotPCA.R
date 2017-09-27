#' Plot Principal Component Analysis (PCA)
#'
#' @rdname plotPCA
#' @name plotPCA
#' @family PCA Utilities
#' @author Michael Steinbaugh
#'
#' @inherit plotTSNE
NULL



# Methods ====
#' @rdname plotPCA
#' @export
setMethod("plotPCA", "seurat", function(
    object,
    interestingGroup = "ident",
    label = TRUE) {
    pca <- fetchPCAData(object)
    .plotPCA(
        pca,
        interestingGroup = interestingGroup,
        label = label)
})
