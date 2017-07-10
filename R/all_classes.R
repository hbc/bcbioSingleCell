#' bcbioSCDataSet
#'
#' [bcbioSCDataSet] is a subclass of
#' [SummarizedExperiment::SummarizedExperiment] designed to store a single-cell
#' RNA-seq analysis. This class contains read counts save as a sparse matrix
#' (`dgCMatrix`), sample barcodes, run metadata, and barcode summary statistics
#' for each sample analyzed.
#'
#' @author Michael Steinbaugh
#'
#' @export
bcbioSCDataSet <- setClass(  # nolint
    "bcbioSCDataSet",
    contains = "SummarizedExperiment",
    representation = representation(callers = "SimpleList"))
setValidity("bcbioSCDataSet", function(object) TRUE)



#' SCSubset
#'
#' [SCSubset] is a subclass of
#' [SummarizedExperiment::SummarizedExperiment] designed to store a subset
#' of single-cell RNA-seq counts filtered by quality control analysis.
#'
#' @author Michael Steinbaugh
#'
#' @export
SCSubset <- setClass("SCSubset", contains = "SummarizedExperiment")  # nolint
