#' bcbioSingleCell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @rdname bcbioSingleCell-package
#' @name bcbioSingleCell-package
#' @docType package
#'
#' @import basejump methods SingleCellExperiment SummarizedExperiment
#'
#' @importClassesFrom Seurat seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @importFrom BiocGenerics plotPCA
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dendsort dendsort
#' @importFrom devtools session_info
#' @importFrom dplyr arrange bind_rows distinct everything filter funs group_by
#'   left_join mutate n pull rename select slice summarize summarize_all ungroup
#' @importFrom ggridges geom_density_ridges
#' @importFrom ggplot2 aes_ aes_string coord_flip element_blank element_rect
#'   element_text expand_limits facet_wrap geom_bar geom_boxplot geom_histogram
#'   geom_hline geom_label geom_line geom_point geom_smooth geom_text
#'   geom_violin geom_vline ggplot ggtitle guide_legend guides labs qplot
#'   scale_color_gradient scale_radius scale_x_log10 scale_x_sqrt
#'   scale_y_continuous scale_y_log10 scale_y_sqrt theme unit xlab xlim ylab
#' @importFrom jsonlite read_json
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom Matrix cBind
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom pbapply pblapply pbsapply
#' @importFrom pheatmap pheatmap
#' @importFrom readr read_csv read_delim read_lines
#' @importFrom rlang .data is_string sym syms
#' @importFrom S4Vectors metadata SimpleList
#' @importFrom scales math_format percent trans_breaks trans_format
#' @importFrom Seurat AddMetaData CreateSeuratObject DarkTheme DotPlot
#'   FeaturePlot FetchData FindVariableGenes JoyPlot NormalizeData ScaleData
#'   VlnPlot
#' @importFrom stats median
#' @importFrom stringr str_detect str_match str_replace str_replace_all
#'   str_split str_subset str_trunc
#' @importFrom tibble column_to_rownames remove_rownames rownames_to_column
#'   tibble
#' @importFrom tidyr gather
#' @importFrom utils globalVariables methods
#' @importFrom viridis inferno scale_color_viridis scale_fill_viridis viridis
NULL

globalVariables(".")
metaPriorityCols <- c("sampleID", "sampleName", "fileName")
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
