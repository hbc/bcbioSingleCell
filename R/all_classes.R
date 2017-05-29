#' Class that contains bcbio run information
#'
#' `bcbioSCDataSet` is a subclass of [SCESet] designed to store a single-cell
#' RNA-seq analysis. This class contains experiment metadata, raw counts,
#' normalilzed counts, and summary statistics for each sample analyzed.
#'
#' Methods for this objects ...
#'
#' `metadata` contains ...
#'
#' @rdname bcbioSCDataSet
#' @keywords internal
#' @aliases bcbioSCDataSet-class
#' @export
bcbioSCDataSet <- setClass("bcbioSCDataSet", contains = "SCESet")
setValidity("bcbioSCDataSet", function(object) TRUE)
