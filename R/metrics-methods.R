#' @name metrics
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit basejump::metrics
#' @inheritParams basejump::params
#'
#' @param recalculate `boolean`. Force recalculation, using primary
#'   [BiocGenerics::counts()] matrix.
#'
#' @examples
#' data(indrops)
#' x <- metrics(indrops, recalculate = TRUE)
#' print(x)
NULL



#' @importFrom basejump metrics
#' @aliases NULL
#' @export
basejump::metrics



metrics.bcbioSingleCell <-  # nolint
    function(object, recalculate = FALSE) {
        validObject(object)
        if (isTRUE(recalculate)) {
            colData <- colData(object)
            metrics <- .metrics.matrix(
                object = counts(object),
                rowRanges = rowRanges(object)
            )
            colData <- colData[
                ,
                setdiff(colnames(colData), colnames(metrics)),
                drop = FALSE
            ]
            colData <- cbind(metrics, colData)
            colData(object) <- colData
            object
        } else {
            metrics(as(object, "SingleCellExperiment"))
        }
    }



#' @rdname metrics
#' @export
setMethod(
    f = "metrics",
    signature = signature("bcbioSingleCell"),
    definition = metrics.bcbioSingleCell
)
