#' bcbioSCDataSet
#'
#' `bcbioSCDataSet` is a subclass of [SummarizedExperiment] designed to store a
#' single-cell RNA-seq analysis. This class contains sparse counts saved as `dgCMatrix`, sample barcodes, other metadata, and summary statistics for each sample analyzed.
#'
#' @rdname bcbioSCDataSet
#' @author Michael Steinbaugh
#'
#' @export
bcbioSCDataSet <- setClass(
    "bcbioSCDataSet",
    contains = "SummarizedExperiment",
    representation = representation(callers = "SimpleList"))
setValidity("bcbioSCDataSet", function(object) TRUE)
