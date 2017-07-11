#' bcbioSinglecell
#'
#' Utility functions for analysis of bcbio-nextgen single-cell RNA-seq data.
#'
#' @import annotables basejump methods SummarizedExperiment S4Vectors
#' @importFrom BiocGenerics counts design
#' @importFrom ggplot2 aes_ coord_flip element_text expand_limits facet_wrap
#'   geom_bar geom_boxplot geom_histogram geom_hline geom_label geom_line
#'   geom_point geom_smooth geom_text geom_violin geom_vline ggplot ggtitle labs
#'   qplot scale_x_log10 scale_x_sqrt scale_y_log10 scale_y_sqrt theme unit xlab
#'   xlim ylab
#' @importFrom Matrix cBind
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom scales math_format trans_breaks trans_format
#' @importFrom Seurat FeaturePlot MultiPlotList Setup SubsetData VlnPlot
#' @importFrom stats median
#' @importFrom utils methods object.size
#' @importClassesFrom Seurat seurat
"_PACKAGE"

globalVariables(".")

bins <- 100L

fail_color <- "red"
pass_color <- "green"
warn_color <- "orange"

meta_priority_cols <- c("sample_id", "sample_name", "file_name")
project_dir_pattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
