#' bcbioSinglecell
#'
#' Utility functions for analysis of bcbio-nextgen single-cell RNA-seq data.
#'
#' @import basejump
#' @import ggplot2
#' @importFrom Matrix colSums readMM writeMM
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom methods as show
#' @importFrom scales math_format trans_breaks trans_format
#' @importFrom stats aggregate cor density median
#' @importFrom utils head
"_PACKAGE"

globalVariables(basejump::globals,
                asNamespace("bcbioSinglecell"),
                add = TRUE)

bins <- 100
fail_color <- "red"
pass_color <- "green"
warn_color <- "orange"
