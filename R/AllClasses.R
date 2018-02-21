setOldClass(Classes = c("grouped_df", "tbl_df", "tibble"))



#' bcbioSingleCell Object Class
#'
#' `bcbioSingleCell` is an extension of `SummarizedExperiment` designed to store
#' bcbio single-cell RNA-seq counts. This class contains read counts save as a
#' sparse matrix (`dgCMatrix`), sample barcodes, run metadata, and barcode
#' summary statistics for each sample analyzed.
#'
#' @note `bcbioSingleCell` extended `SummarizedExperiment` until v0.0.31 of
#' this package, when we migrated over to using `SingleCellExperiment`.
#'
#' @author Michael Steinbaugh
#'
#' @slot bcbio `SimpleList` containing additional bcbio run data with dimensions
#' that don't match the count matrix. This is currently used to store all
#' unfiltered cellular barcodes for quality control analysis.
#'
#' @seealso
#' - [SummarizedExperiment::SummarizedExperiment()].
#' - `.S4methods(class = "bcbioSingleCell")`.
#'
#' @export
bcbioSingleCell <- setClass(
    "bcbioSingleCell",
    contains = "SummarizedExperiment",
    slots = c(bcbio = "SimpleList")
)
setValidity(
    "bcbioSingleCell",
    function(object) {
        TRUE
    }
)



# Legacy classes ===============================================================
#' `bcbioSCDataSet`
#'
#' This class will be deprecated in favor of [bcbioSingleCell] in a future
#' release.
#'
#' @author Michael Steinbaugh
#' @keywords internal
#'
#' @slot callers [SimpleList] containing additional bcbio run data with
#'   dimensions that don't match the count matrix. This is currently used to
#'   store all unfiltered cellular barcodes for quality control analysis.
#'
#' @export
bcbioSCDataSet <- setClass(
    "bcbioSCDataSet",
    contains = "SummarizedExperiment",
    slots = c(callers = "SimpleList"))
setValidity("bcbioSCDataSet", function(object) TRUE)



#' `bcbioSCFiltered`
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
