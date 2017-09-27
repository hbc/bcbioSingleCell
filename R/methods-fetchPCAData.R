#' Fetch PCA Dimensions and Cellular Metadata
#'
#' @rdname fetchPCAData
#' @name fetchPCAData
#' @family PCA Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return [data.frame] of PCA points and metadata for each cell.
#'
#' @seealso [Seurat::PCAPlot()].
#'
#' @examples
#' \dontrun{
#' data(seurat)
#' pca <- fetchPCAData(seurat)
#' }
NULL



# Methods ====
#' @rdname fetchPCAData
#' @export
setMethod("fetchPCAData", "seurat", function(object) {
    .fetchDimDataSeurat(
        object,
        dimCode = c(x = "PC1", y = "PC2")
    )
})
