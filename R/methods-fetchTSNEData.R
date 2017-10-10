#' Fetch t-SNE Locations and Cellular Metadata
#'
#' @rdname fetchTSNEData
#' @name fetchTSNEData
#' @family t-SNE Utilities
#' @author Rory Kirchner, Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @return [data.frame] of t-SNE points and metadata for each cell.
#'
#' @seealso [Seurat::TSNEPlot()].
#'
#' @examples
#' \dontrun{
#' data(seurat)
#' tsne <- fetchTSNEData(seurat)
#' }
NULL



# Methods ====
#' @rdname fetchTSNEData
#' @export
setMethod(
    "fetchTSNEData",
    signature("seurat"),
    function(object) {
        .fetchDimDataSeurat(
            object,
            dimCode = c(x = "tSNE_1", y = "tSNE_2")
        )
    })
