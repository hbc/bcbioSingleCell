#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom BiocGenerics counts updateObject
#' @importFrom BiocParallel bplapply bpmapply
#' @importFrom IRanges DataFrameList
#' @importFrom S4Vectors DataFrame SimpleList cbind do.call lapply metadata
#'   metadata<-
#' @importFrom SummarizedExperiment assayNames assay assays assays<- colData
#'   colData<- rowData rowData<- rowRanges rowRanges<-
#' @importFrom acidplots acid_geom_abline acid_geom_label
#'   acid_geom_label_average acid_geom_label_repel
#' @importFrom basejump as_tibble calculateMetrics camelCase cell2sample
#'   detectLanes droplevels emptyRanges formalsList import importSampleData
#'   interestingGroups interestingGroups<- leftJoin makeDimnames
#'   makeGRangesFromEnsembl makeGRangesFromGFF makeLabel makeNames
#'   makeSingleCellExperiment mapCellsToSamples markdownPlots matchArgsToDoCall
#'   matchInterestingGroups metrics metricsCols minimalSampleData printString
#'   realpath sampleData sampleNames separator showSlotInfo standardizeCall
#' @importFrom bcbioBase getBarcodeCutoffFromCommands getGTFFileFromYAML
#'   getLevelFromCommands getSampleDataFromYAML getUMITypeFromCommands
#'   importDataVersions importProgramVersions projectDir runDate sampleDirs
#' @importFrom cli cat_line cli_alert cli_alert_success cli_alert_warning cli_h1
#'   cli_h2 cli_text
#' @importFrom ggplot2 aes facet_wrap geom_boxplot geom_histogram geom_step
#'   geom_violin ggplot labs scale_x_continuous scale_y_continuous stat_ecdf
#' @importFrom ggridges geom_density_ridges
#' @importFrom goalie allAreDirectories allAreFiles areDisjointSets areSetEqual
#'   assert hasLength hasNames hasRownames hasValidDimnames isADirectory isAFile
#'   isAURL isAny isCharacter isDirectory isFile isFlag isGGScale isInt isString
#'   isSubset validate validateClasses
#' @importFrom graphics hist
#' @importFrom methods as as<- is new setClass show slot slot<- validObject
#'   .hasSlot
#' @importFrom rlang !! sym syms
#' @importFrom stringr str_extract
#' @importFrom utils capture.output globalVariables packageVersion
NULL
