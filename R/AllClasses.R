# bcbioSingleCell ==============================================================
#' `bcbioSingleCell` Class
#'
#' `bcbioSingleCell` is an S4 class that extends `SingleCellExperiment`, and is
#' designed to store a bcbio single-cell RNA-seq analysis. This class contains
#' read counts saved as a sparse matrix (`dgCMatrix`), sample metadata, and cell
#' quality control metrics.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @seealso [bcbioSingleCell()].
setClass(
    Class = "bcbioSingleCell",
    contains = "SingleCellExperiment"
)



# CellRanger ===================================================================
#' `CellRanger` Class
#'
#' Extends `SingleCellExperiment`, with additional validity checks on the
#' `metadata()` slot.
#'
#' @family S4 Classes
#' @author Michael Steinbaugh
#' @export
#'
#' @seealso [CellRanger()].
setClass(
    Class = "CellRanger",
    contains = "SingleCellExperiment"
)
