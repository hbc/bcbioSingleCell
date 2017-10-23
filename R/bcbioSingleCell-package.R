#' bcbioSingleCell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @rdname bcbioSingleCell-package
#' @name bcbioSingleCell-package
#' @docType package
#'
#' @import basejump methods SummarizedExperiment
#'
#' @importClassesFrom Seurat seurat
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom utils globalVariables
#'
#' @importFrom dplyr arrange bind_rows distinct everything filter funs group_by
#'   left_join mutate n pull rename select slice summarize summarize_all ungroup
#' @importFrom ggplot2 aes_ aes_string coord_flip element_blank element_line
#'   element_rect element_text expand_limits facet_wrap geom_bar geom_boxplot
#'   geom_histogram geom_hline geom_label geom_line geom_point geom_smooth
#'   geom_text geom_violin geom_vline ggplot ggtitle guide_colorbar
#'   guide_legend guides labs qplot scale_color_gradient scale_radius
#'   scale_x_log10 scale_x_sqrt scale_y_continuous scale_y_log10 scale_y_sqrt
#'   theme unit xlab xlim ylab
#' @importFrom stringr str_detect str_match str_replace str_replace_all
#'   str_split str_subset str_trunc
NULL

utils::globalVariables(".")
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
