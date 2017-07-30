#' bcbioSinglecell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @import basejump methods SummarizedExperiment S4Vectors
#' @importFrom BiocGenerics counts design
#' @importFrom ggplot2 aes_ coord_flip element_blank element_text expand_limits
#'   facet_wrap geom_bar geom_boxplot geom_histogram geom_hline geom_label
#'   geom_line geom_point geom_smooth geom_text geom_violin geom_vline ggplot
#'   ggtitle labs qplot scale_x_log10 scale_x_sqrt scale_y_log10 scale_y_sqrt
#'   theme unit xlab xlim ylab
#' @importFrom Matrix cBind
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom scales math_format trans_breaks trans_format
#' @importFrom Seurat AddMetaData CreateSeuratObject FeaturePlot NormalizeData
#'   VlnPlot
#' @importFrom stats median
#' @importFrom utils methods object.size
#' @importClassesFrom Seurat seurat
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
