#' bcbioSinglecell
#'
#' Utility functions for analysis of bcbio-nextgen single-cell RNA-seq data.
#'
#' @keywords internal
#'
#' @import basejump
#' @import Biobase
#' @import ggplot2
#' @import methods
#' @import SummarizedExperiment
#' @importFrom Matrix colSums readMM rowSums writeMM
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom scales math_format trans_breaks trans_format
#' @importFrom stats aggregate median
#' @importFrom utils methods
"_PACKAGE"

globalVariables(".")

bins <- 100
fail_color <- "red"
pass_color <- "green"
warn_color <- "orange"
