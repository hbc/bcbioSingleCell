# TODO Add validity checks



setOldClass(Classes = c("grouped_df", "tbl_df", "tibble"))



#' bcbioSingleCell Class
#'
#' `bcbioSingleCell` extends `SingleCellExperiment` and designed to store a
#' bcbio single-cell RNA-seq analysis. This class contains read counts saved as
#' a sparse matrix (`dgCMatrix`), sample metadata, and cell quality control
#' metrics.
#'
#' @note `bcbioSingleCell` extended `SummarizedExperiment` prior to v0.1.0,
#'   where we migrated to `SingleCellExperiment`.
#'
#' @author Michael Steinbaugh
#'
#' @slot bcbio `SimpleList` containing additional bcbio run data with dimensions
#' that don't match the count matrix. This is currently used to store all
#' unfiltered cellular barcodes for quality control analysis.
#'
#' @seealso
#' - [loadSingleCell()], [loadCellRanger()].
#' - [SingleCellExperiment::SingleCellExperiment()].
#' - `.S4methods(class = "bcbioSingleCell")`.
#'
#' @export
bcbioSingleCell <- setClass(
    "bcbioSingleCell",
    contains = "SingleCellExperiment"
)



# Validity =====================================================================
setValidity(
    "bcbioSingleCell",
    function(object) {
        TRUE
    }
)
# object@bcbio$cellularBarcodes needs to go in metadata
