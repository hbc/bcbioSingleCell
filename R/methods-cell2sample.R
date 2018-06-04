#' Cell to Sample Mappings
#'
#' @name cell2sample
#' @family Data Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return Named `factor` containing sample IDs as the levels and cellular
#'   barcode IDs as the names.
#'
#' @examples
#' # SingleCellExperiment ====
#' x <- cell2sample(cellranger_small)
#' table(x)
NULL



# Methods ======================================================================
#' @rdname cell2sample
#' @export
setMethod(
    "cell2sample",
    signature("SingleCellExperiment"),
    function(object) {
        validObject(object)
        stash <- metadata(object)[["cell2sample"]]
        if (
            is.factor(stash) &&
            identical(colnames(object), names(stash))
        ) {
            return(stash)
        } else {
            mapCellsToSamples(
                cells = colnames(object),
                samples = rownames(sampleData(object))
            )
        }
    }
)



#' @rdname seurat-SingleCellExperiment
#' @export
setMethod(
    "cell2sample",
    signature("seurat"),
    getMethod("cell2sample", "SingleCellExperiment")
)
