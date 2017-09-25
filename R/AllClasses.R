setOldClass(c("grouped_df", "tbl_df", "tibble"))



#' bcbioSCDataSet
#'
#' [bcbioSCDataSet] is an extension of [SummarizedExperiment] designed to store
#' a single-cell RNA-seq analysis. This class contains read counts save as a
#' sparse matrix (`dgCMatrix`), sample barcodes, run metadata, and barcode
#' summary statistics for each sample analyzed.
#'
#' @author Michael Steinbaugh
#'
#' @slot bcbio [SimpleList] containing additional bcbio run data with dimensions
#' that don't match the count matrix. This is currently used to store all
#' unfiltered cellular barcodes for quality control analysis.
#'
#' @export
bcbioSCDataSet <- setClass(
    "bcbioSCDataSet",
    contains = "SummarizedExperiment",
    slots = c(bcbio = "SimpleList"))
setValidity("bcbioSCDataSet", function(object) TRUE)



#' bcbioSCFiltered
#'
#' [bcbioSCFiltered] is a subclass of [SummarizedExperiment] designed to store
#' single-cell RNA-seq counts of cells filtered by quality control analysis.
#'
#' @author Michael Steinbaugh
#'
#' @export
bcbioSCFiltered <- setClass(
    "bcbioSCFiltered",
    contains = "SummarizedExperiment")
setValidity("bcbioSCFiltered", function(object) TRUE)
