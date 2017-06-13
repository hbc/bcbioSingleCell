#' Class that contains bcbio run information
#'
#' `bcbioSCDataSet` is a subclass of [SummarizedExperiment] designed to store a
#' single-cell RNA-seq analysis. This class contains experiment metadata, raw
#' counts, normalilzed counts, and summary statistics for each sample analyzed.
#'
#' Methods for this objects ...
#'
#' `metadata` contains ...
#'
#' @rdname bcbioSCDataSet
#' @keywords internal
#' @aliases bcbioSCDataSet-class
#' @export
bcbioSCDataSet <- setClass(
    "bcbioSCDataSet", contains = "SummarizedExperiment")
# setValidity("bcbioSCDataSet", function(object) TRUE)



# Constructor
.bcbioRnaDataSet <- function(se, run) {
    if (!is(se, "SummarizedExperiment")) {
        if (is(se, "SummarizedExperiment0")) {
            se <- as(se, "SummarizedExperiment")
        } else if (is(se, "SummarizedExperiment")) {
            # Only to help transition from SummarizedExperiment to new
            # RangedSummarizedExperiment objects, remove once transition is
            # complete
            se <- as(se, "SummarizedExperiment")
        } else {
            stop("'se' must be a SummarizedExperiment object")
        }
    }
    bcb <- new("bcbioRnaDataSet", se)
    bcb
}
