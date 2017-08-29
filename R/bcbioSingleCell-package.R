#' bcbioSingleCell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @import methods
#' @importClassesFrom Biobase AnnotatedDataFrame
#' @importClassesFrom monocle CellDataSet
#' @importClassesFrom Seurat seurat
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#' @importFrom basejump camel detectOrganism gene2symbolFromGTF mdHeader mdList
#'   prepareSE prepareTemplate readFileByExtension readGTF revcomp tidy_filter
#'   tidy_select tx2geneFromGTF
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dplyr arrange bind_rows distinct everything group_by left_join
#'   mutate n pull rename summarize top_n
#' @importFrom ggplot2 aes_ coord_flip element_blank element_text expand_limits
#'   facet_wrap geom_bar geom_boxplot geom_histogram geom_hline geom_label
#'   geom_line geom_point geom_smooth geom_text geom_violin geom_vline ggplot
#'   ggtitle labs qplot scale_x_log10 scale_x_sqrt scale_y_continuous
#'   scale_y_log10 scale_y_sqrt theme unit xlab xlim ylab
#' @importFrom jsonlite read_json
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_rownames
#' @importFrom Matrix cBind
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom monocle newCellDataSet
#' @importFrom pbapply pblapply pbsapply
#' @importFrom readr read_csv read_delim read_lines
#' @importFrom rlang .data is_string sym syms
#' @importFrom S4Vectors metadata SimpleList
#' @importFrom scales math_format percent trans_breaks trans_format
#' @importFrom Seurat AddMetaData CreateSeuratObject FeaturePlot
#'   FindVariableGenes NormalizeData ScaleData
#'   VlnPlot
#' @importFrom stats median
#' @importFrom stringr str_detect str_match str_replace str_replace_all
#'   str_split str_subset
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom tibble column_to_rownames remove_rownames rownames_to_column
#'   tibble
#' @importFrom utils globalVariables methods
#' @importFrom viridis scale_color_viridis scale_fill_viridis
"_PACKAGE"

globalVariables(".")

bins <- 100L

# Quality control plot colors
qcCutoffColor <- "black"
qcFailColor <- "red"
qcPassColor <- "green"
qcWarnColor <- "orange"

# Quality control label appearance
qcLabelAlpha <- 0.75
qcLabelColor <- "white"
qcLabelFill <- "black"
qcLabelFontface <- "bold"
qcLabelPadding <- unit(0.2, "lines")
qcLabelSize <- NA

# Quality control line appearance
qcLineAlpha <- 0.75
qcLineSize <- 1.5
qcLineType <- "longdash"

# Plot label separator
labelSep <- ": "

metaPriorityCols <- c("sampleID", "sampleName", "fileName")
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"

url <- "http://bioinformatics.sph.harvard.edu/bcbioSingleCell"
