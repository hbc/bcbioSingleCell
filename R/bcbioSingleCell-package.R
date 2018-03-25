#' bcbioSingleCell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @rdname bcbioSingleCell-package
#' @name bcbioSingleCell-package
#' @docType package
#'
#' @import methods S4Vectors
#'
#' @importClassesFrom Seurat seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @importFrom Biobase rowMedians
#' @importFrom GenomicFeatures makeTxDbFromGFF
#' @importFrom Matrix cBind readMM
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom SingleCellExperiment isSpike isSpike<- SingleCellExperiment
#' @importFrom Seurat CreateSeuratObject FetchData
#' @importFrom SummarizedExperiment assay assays colData rowData rowRanges
#' @importFrom basejump assignAndSaveData camel convertUCSCBuildToEnsembl
#'   detectOrganism dynamicPlotlist ensembl initializeDirectory makeNames
#'   markdownHeader midnightTheme readFileByExtension readYAML
#'   sanitizeSampleData stripTranscriptVersions
#' @importFrom bcbioBase prepareSummarizedExperiment readDataVersions
#'   readLogFile readProgramVersions readSampleMetadataFile sampleYAMLMetadata
#'   uniteInterestingGroups
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dplyr arrange bind_rows group_by group_vars left_join mutate
#'   mutate_all mutate_if n pull select_if summarize summarize_all ungroup
#' @importFrom edgeR calcNormFactors DGEList glmFit
#' @importFrom ggplot2 aes_ aes_string coord_flip element_blank element_line
#'   element_rect element_text expand_limits facet_wrap geom_bar geom_boxplot
#'   geom_histogram geom_hline geom_label geom_line geom_point geom_smooth
#'   geom_text geom_violin geom_vline ggplot ggtitle guide_colorbar guide_legend
#'   guides labs qplot scale_radius scale_x_log10 scale_x_sqrt
#'   scale_y_continuous scale_y_log10 scale_y_sqrt theme xlab xlim ylab
#' @importFrom ggridges geom_density_ridges
#' @importFrom graphics hist
#' @importFrom grid unit
#' @importFrom jsonlite read_json
#' @importFrom magrittr %>% set_colnames set_names set_rownames
#' @importFrom parallel mclapply mcmapply
#' @importFrom pbapply pblapply
#' @importFrom readr read_lines read_tsv
#' @importFrom rlang !! !!! .data abort inform sym syms warn
#' @importFrom scales percent
#' @importFrom stats model.matrix reorder
#' @importFrom stringr str_extract str_match str_pad str_split
#' @importFrom tibble as_tibble column_to_rownames remove_rownames
#'   rownames_to_column tibble
#' @importFrom tidyr gather
#' @importFrom utils globalVariables packageVersion
#' @importFrom zingeR glmWeightedF zeroWeightsLS
NULL



globalVariables(".")
packageVersion <- packageVersion("bcbioSingleCell")

# Trailing number is to match cellranger output
barcodePattern <- ")_([ACGT_]{6,})(_[0-9]+)?$"
lanePattern <- "_L(\\d{3})"
metadataPriorityCols <- c("sampleID", "sampleName", "description")
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"

sepBar <- "============================================================"

validMedianGeom <- c(
    "boxplot",
    "ridgeline",
    "violin"
)
validQCGeomFlip <- c(
    "boxplot",
    "violin"
)
