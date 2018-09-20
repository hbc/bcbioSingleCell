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
#' x <- cell2sample(indrops_small)
#' table(x)
NULL



.cell2sample.SCE <-  # nolint
    function(object) {
        validObject(object)
        stash <- metadata(object)[["cell2sample"]]
        if (
            is.factor(stash) &&
            identical(colnames(object), names(stash))
        ) {
            return(stash)
        }
        cells <- colnames(object)
        samples <- rownames(sampleData(object))
        if (is.null(samples)) {
            samples <- "unknown"
        }
        .mapCellsToSamples(cells = cells, samples = samples)
    }



#' @rdname cell2sample
#' @export
setMethod(
    f = "cell2sample",
    signature = signature("SingleCellExperiment"),
    definition = .cell2sample.SCE
)
