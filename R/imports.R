#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom BiocGenerics cbind counts do.call updateObject
#' @importFrom BiocParallel bplapply bpmapply
#' @importFrom IRanges DataFrameList
#' @importFrom Matrix readMM
#' @importFrom S4Vectors DataFrame SimpleList metadata metadata<-
#' @importFrom SummarizedExperiment assayNames assay assays assays<- colData
#'   colData<- rowData rowData<- rowRanges rowRanges<-
#' @importFrom acidplots acid_geom_abline acid_geom_label
#'   acid_geom_label_average acid_geom_label_repel
#' @importFrom basejump as_tibble calculateMetrics camelCase cell2sample
#'   detectLanes droplevels emptyRanges formalsList import interestingGroups
#'   interestingGroups<- left_join makeDimnames makeGRangesFromEnsembl
#'   makeGRangesFromGFF makeLabel makeNames makeSingleCellExperiment
#'   mapCellsToSamples markdownPlots matchArgsToDoCall matchInterestingGroups
#'   metrics minimalSampleData printString realpath sampleData
#'   sampleNames separator showSlotInfo
#' @importFrom bcbioBase getBarcodeCutoffFromCommands getGTFFileFromYAML
#'   getLevelFromCommands getSampleDataFromYAML getUMITypeFromCommands
#'   projectDir readDataVersions readProgramVersions readSampleData runDate
#'   sampleDirs
#' @importFrom ggplot2 aes facet_wrap geom_boxplot geom_histogram geom_step
#'   geom_violin ggplot labs scale_x_continuous scale_y_continuous stat_ecdf
#' @importFrom ggridges geom_density_ridges
#' @importFrom goalie allAreDirectories allAreFiles areDisjointSets areSetEqual
#'   assert hasLength hasNames hasRownames hasValidDimnames isADirectory isAFile
#'   isAURL isAny isCharacter isDirectory isFile isFlag isGGScale isInt
#'   isNonEmpty isString isSubset validate validateClasses
#' @importFrom graphics hist
#' @importFrom methods as as<- is new setClass show slot slot<- validObject
#'   .hasSlot
#' @importFrom rlang !! sym syms
#' @importFrom stringr str_extract
#' @importFrom tibble as_tibble tibble
#' @importFrom utils capture.output globalVariables packageVersion
#'
#'
#'
#' @importFrom dplyr bind_rows group_by mutate mutate_all mutate_if
#' @importFrom magrittr %>%
NULL
