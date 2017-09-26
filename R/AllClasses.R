#' bcbioSingleCell
#'
#' [bcbioSingleCell] is an extension of [SingleCellExperiment] designed to store
#' bcbio single-cell RNA-seq counts. This class contains read counts save as a
#' sparse matrix (`dgCMatrix`), sample barcodes, run metadata, and barcode
#' summary statistics for each sample analyzed.
#'
#' @author Michael Steinbaugh
#'
#' @slot bcbio [SimpleList] containing additional bcbio run data with dimensions
#' that don't match the count matrix. This is currently used to store all
#' unfiltered cellular barcodes for quality control analysis.
#'
#' @seealso `.S4methods(class = "bcbioSingleCell")`
#'
#' @export
bcbioSingleCell <- setClass(
    "bcbioSingleCell",
    contains = "SingleCellExperiment",
    slots = c(bcbio = "SimpleList"))
setValidity("bcbioSingleCell", function(object) TRUE)



# Future Class Deprecations ====
#' bcbioSCDataSet
#'
#' This class will be deprecated in favor of [bcbioSingleCell] in a future
#' release.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @slot callers [SimpleList] containing additional bcbio run data with dimensions
#' that don't match the count matrix. This is currently used to store all
#' unfiltered cellular barcodes for quality control analysis.
#'
#' @export
bcbioSCDataSet <- setClass(
    "bcbioSCDataSet",
    contains = "SummarizedExperiment",
    slots = c(callers = "SimpleList"))
setValidity("bcbioSCDataSet", function(object) TRUE)



#' bcbioSCFiltered
#'
#' This class will be deprecated in favor of [bcbioSingleCell] in a future
#' release.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @export
bcbioSCFiltered <- setClass(
    "bcbioSCFiltered",
    contains = "SummarizedExperiment")
setValidity("bcbioSCFiltered", function(object) TRUE)



# Miscellaneous Class Definitions ====
# Need to register tibble classes for S4 method support
setOldClass(
    Classes = c("grouped_df", "tbl_df", "tibble"))
# Remove this once the old `bcbioSC*` classes are formally deprecated
setClassUnion(
    name = "bcbioSingleCellANY",
    members = c("bcbioSingleCell",
                "bcbioSCDataSet",
                "bcbioSCFiltered"))
setClassUnion(
    name = "bcbioSingleCellLegacy",
    members = c("bcbioSCDataSet",
                "bcbioSCFiltered"))
