#' Fetch PCA Dimensions and Cellular Metadata
#'
#' @name fetchPCAData
#' @family PCA Utilities
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `data.frame` of PCA points and metadata for each cell.
#'
#' @seealso [Seurat::PCAPlot()].
#'
#' @examples
#' # seurat ====
#' fetchPCAData(pbmc_small) %>% glimpse()
NULL



# Methods ======================================================================
#' @rdname fetchPCAData
#' @export
setMethod(
    "fetchPCAData",
    signature("seurat"),
    function(object) {
        .fetchDRData.seurat(
            object,
            dimCode = c(x = "PC1", y = "PC2")
        )
    })
