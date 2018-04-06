#' bcbioSingleCell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @name bcbioSingleCell-package
#' @keywords internal
#'
#' @importClassesFrom Seurat seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Biobase rowMedians
#' @importFrom GenomicFeatures genes makeTxDbFromGFF transcripts
#' @importFrom Matrix colSums readMM rowMeans rowSums
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom S4Vectors aggregate complete.cases metadata metadata<- na.omit
#' @importFrom SingleCellExperiment isSpike isSpike<- SingleCellExperiment
#' @importFrom Seurat CreateSeuratObject FetchData
#' @importFrom SummarizedExperiment assay assays colData rowData rowRanges
#' @importFrom basejump assignAndSaveData camel convertUCSCBuildToEnsembl
#'   detectOrganism dynamicPlotlist initializeDirectory makeGRangesFromEnsembl
#'   makeGRangesFromGFF makeNames makeTx2geneFromGFF markdownHeader
#'   midnightTheme readFileByExtension readYAML sanitizeSampleData
#'   tx2geneFromGFF
#' @importFrom bcbioBase flatFiles prepareSummarizedExperiment readDataVersions
#'   readLog readProgramVersions readSampleData readTx2gene
#'   sampleYAMLMetadata uniteInterestingGroups
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dplyr arrange bind_rows group_by group_vars left_join matches
#'   mutate mutate_all mutate_if n select_if summarize summarize_all ungroup
#' @importFrom edgeR calcNormFactors DGEList estimateDisp glmFit
#' @importFrom ggplot2 aes_ aes_string coord_flip element_blank element_line
#'   element_rect element_text expand_limits facet_wrap geom_bar geom_boxplot
#'   geom_histogram geom_hline geom_label geom_line geom_point geom_smooth
#'   geom_text geom_violin geom_vline ggtitle guide_colorbar guide_legend guides
#'   labs qplot scale_radius scale_x_continuous scale_y_continuous stat_ecdf
#'   theme xlab xlim ylab
#' @importFrom ggridges geom_density_ridges
#' @importFrom graphics hist
#' @importFrom grid unit
#' @importFrom jsonlite read_json
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_names set_rownames
#' @importFrom methods .hasSlot as is new show slot slot<- validObject
#' @importFrom parallel mclapply mcmapply
#' @importFrom pbapply pblapply
#' @importFrom readr read_lines read_tsv
#' @importFrom rlang !! !!! .data abort inform sym syms warn
#' @importFrom scales percent
#' @importFrom stats as.formula fitted median model.matrix predict relevel
#'   reorder smooth.spline
#' @importFrom stringr str_extract str_match str_pad str_split
#' @importFrom tibble as_tibble column_to_rownames remove_rownames
#'   rownames_to_column tibble
#' @importFrom tidyr gather
#' @importFrom utils globalVariables packageVersion
#' @importFrom zingeR glmWeightedF zeroWeightsLS
NULL
