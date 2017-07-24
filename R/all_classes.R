#' bcbioSCDataSet
#'
#' [bcbioSCDataSet] is a subclass of [SummarizedExperiment] designed to store a
#' single-cell RNA-seq analysis. This class contains read counts save as a
#' sparse matrix (`dgCMatrix`), sample barcodes, run metadata, and barcode
#' summary statistics for each sample analyzed.
#'
#' @author Michael Steinbaugh
#'
#' @export
bcbioSCDataSet <- setClass(  # nolint
    "bcbioSCDataSet",
    contains = "SummarizedExperiment",
    slots = c(callers = "SimpleList"))
setValidity("bcbioSCDataSet", function(object) TRUE)



#' bcbioSCSubset
#'
#' [bcbioSCSubset] is a subclass of [SummarizedExperiment] designed to store a
#' subset of single-cell RNA-seq counts filtered by quality control analysis.
#'
#' @author Michael Steinbaugh
#'
#' @export
bcbioSCSubset <- setClass(  # nolint
    "bcbioSCSubset",
    contains = "SummarizedExperiment")
setValidity("bcbioSCSubset", function(object) TRUE)
