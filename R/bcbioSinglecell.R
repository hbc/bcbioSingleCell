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
#' @import S4Vectors
#' @importClassesFrom scater SCESet
#' @importFrom Matrix readMM writeMM
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom scales math_format trans_breaks trans_format
#' @importFrom stats density median
"_PACKAGE"

globalVariables(basejump::globals, asNamespace("bcbioSinglecell"), add = TRUE)

bins <- 100
fail_color <- "red"
pass_color <- "green"
warn_color <- "orange"
