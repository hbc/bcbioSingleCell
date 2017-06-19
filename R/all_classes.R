#' bcbioSCDataSet
#'
#' `bcbioSCDataSet` is a subclass of [SummarizedExperiment] designed to store a
#' single-cell RNA-seq analysis. This class contains experiment metadata, raw
#' counts, and summary statistics for each sample analyzed.
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
    "bcbioSCDataSet",
    contains = "SummarizedExperiment",
    representation = representation(
        callers = "SimpleList"))
setValidity("bcbioSCDataSet", function(object) TRUE)
