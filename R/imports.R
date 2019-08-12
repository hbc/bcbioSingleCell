#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom BiocGenerics cbind counts do.call updateObject
#' @importFrom BiocParallel bplapply bpmapply
#' @importFrom IRanges DataFrameList
#' @importFrom Matrix colSums readMM rowSums
#' @importFrom S4Vectors SimpleList as.data.frame as.matrix mcols
#'   mcols<- merge metadata metadata<-
#' @importFrom SingleCellExperiment isSpike
#' @importFrom SummarizedExperiment assayNames assays colData colData<- rowData
#'   rowData<- rowRanges rowRanges<-
#' @importFrom acidplots acid_geom_bar acid_scale_y_continuous_nopad synesthesia
#' @importFrom basejump as_tibble calculateMetrics camelCase cell2sample
#'   detectLanes emptyRanges formalsList import interestingGroups
#'   interestingGroups<- left_join makeDimnames makeGRangesFromEnsembl
#'   makeGRangesFromGFF makeLabel makeNames makeSingleCellExperiment makeTitle
#'   mapCellsToSamples markdownHeader markdownPlots matchArgsToDoCall
#'   matchInterestingGroups metrics minimalSampleData prepareTemplate
#'   printString relevel realpath sampleData sampleNames separator showSlotInfo
#' @importFrom bcbioBase getBarcodeCutoffFromCommands getGTFFileFromYAML
#'   getLevelFromCommands getSampleDataFromYAML getUMITypeFromCommands
#'   projectDir readDataVersions readProgramVersions readSampleData runDate
#'   sampleDirs
#' @importFrom cowplot plot_grid
#' @importFrom dplyr arrange bind_rows desc filter group_by left_join mutate
#'   mutate_all mutate_if n pull slice summarise
#' @importFrom ggplot2 aes facet_wrap geom_bar geom_boxplot geom_histogram
#'   geom_hline geom_line geom_point geom_smooth geom_step geom_violin ggplot
#'   labs scale_x_continuous scale_y_continuous stat_ecdf theme
#' @importFrom ggridges geom_density_ridges
#' @importFrom goalie allAreDirectories allAreFiles allAreHexColors
#'   areDisjointSets areSetEqual assert hasLength hasNames hasRows hasRownames
#'   hasValidDimnames isADirectory isAFile isAURL isAny isCharacter isDirectory
#'   isFile isFlag isGGScale isHeaderLevel isInLeftOpenRange isInRange
#'   isInRightOpenRange isInt isNonEmpty isNonNegative isNumber isPositive
#'   isString isSubset validate validateClasses
#' @importFrom graphics hist
#' @importFrom magrittr %<>% %>% set_names
#' @importFrom methods as as<- is new setClass show slot slot<- validObject
#'   .hasSlot
#' @importFrom acidplots acid_geom_abline acid_geom_label
#'   acid_geom_label_average acid_geom_label_repel plotQC plotZerosVsDepth
#' @importFrom readr read_lines read_tsv
#' @importFrom rlang !! := UQ sym syms
#' @importFrom scales percent
#' @importFrom stats ecdf
#' @importFrom stringr str_extract str_pad
#' @importFrom tibble as_tibble tibble
#' @importFrom utils capture.output globalVariables packageVersion
NULL
