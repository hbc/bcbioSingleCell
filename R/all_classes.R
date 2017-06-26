#' bcbioSCDataSet
#'
#' `bcbioSCDataSet` is a subclass of [SummarizedExperiment] designed to store a
#' single-cell RNA-seq analysis. This class contains read counts save as a
#' sparse matrix (`dgCMatrix`), sample barcodes, run metadata, and barcode
#' summary statistics for each sample analyzed.
#'
#' @rdname bcbioSCDataSet
#' @author Michael Steinbaugh
#'
#' @export
bcbioSCDataSet <- setClass(  # nolint
    "bcbioSCDataSet",
    contains = "SummarizedExperiment",
    representation = representation(callers = "SimpleList"))
setValidity("bcbioSCDataSet", function(object) TRUE)



#' bcbioSCSubset
#'
#' `bcbioSCSubset` is a subclass of [ExpressionSet].
#'
#' @rdname bcbioSCSubset
#' @author Michael Steinbaugh
#'
#' @export
bcbioSCSubset <- setClass(
    "bcbioSCSubset",
    contains = "ExpressionSet",
    representation = representation(callers = "SimpleList"))
setValidity("bcbioSCSubset", function(object) TRUE)
