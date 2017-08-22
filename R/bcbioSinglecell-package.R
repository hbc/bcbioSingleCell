#' bcbioSinglecell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @importClassesFrom Seurat seurat
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom basejump camel gene2symbolFromGTF mdHeader mdList packageSE
#'   readFileByExtension readGTF revcomp tidy_filter tidy_select tx2geneFromGTF
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dplyr arrange bind_rows distinct everything group_by left_join
#'   mutate n pull rename summarize top_n
#' @importFrom ggplot2 aes_ coord_flip element_blank element_text expand_limits
#'   facet_wrap geom_bar geom_boxplot geom_histogram geom_hline geom_label
#'   geom_line geom_point geom_smooth geom_text geom_violin geom_vline ggplot
#'   ggtitle labs qplot scale_x_log10 scale_x_sqrt scale_y_log10 scale_y_sqrt
#'   theme unit xlab xlim ylab
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom Matrix cBind
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom pbapply pblapply
#' @importFrom readr read_csv read_delim read_lines
#' @importFrom rlang .data is_string sym syms
#' @importFrom S4Vectors metadata SimpleList
#' @importFrom scales math_format trans_breaks trans_format
#' @importFrom Seurat AddMetaData CreateSeuratObject FeaturePlot NormalizeData
#'   VlnPlot
#' @importFrom stats median
#' @importFrom stringr str_detect str_match str_replace str_replace_all
#'   str_split str_subset
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble column_to_rownames remove_rownames rownames_to_column
#'   tibble
#' @importFrom utils methods
"_PACKAGE"

globalVariables(".")

bins <- 100L

# Quality control plot colors
qcFailColor <- "red"
qcPassColor <- "green"
qcWarnColor <- "orange"
qcLabelAlpha <- 0.5
qcLineAlpha <- 0.5
qcLineSize <- 2L

# Plot label separator
labelSep <- ": "

metaPriorityCols <- c("sampleID", "sampleName", "fileName")
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
