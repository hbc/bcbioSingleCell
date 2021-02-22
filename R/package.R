#' bcbioSingleCell
#'
#' Import and analyze [bcbio](https://bcbio-nextgen.readthedocs.io/) single-cell
#' RNA-seq data.
#'
#' @aliases NULL
#' @keywords internal
#'
#' @importClassesFrom basejump SingleCellExperiment
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom AcidPlots !! acid_geom_abline acid_geom_label
#'   acid_geom_label_average acid_geom_label_repel sym syms
#' @importFrom BiocParallel bplapply bpmapply bpparam
#' @importFrom basejump DataFrame DataFrameList SimpleList assayNames assay
#'   assays assays<- as_tibble calculateMetrics camelCase capture.output cbind
#'   cell2sample colData colData<- counts detectLanes do.call droplevels
#'   emptyRanges formalsList import importSampleData interestingGroups
#'   interestingGroups<- lapply leftJoin makeDimnames makeGRangesFromEnsembl
#'   makeGRangesFromGFF makeLabel makeNames makeSingleCellExperiment
#'   mapCellsToSamples markdownPlots matchInterestingGroups metadata metadata<-
#'   metrics metricsCols minimalSampleData packageName packageVersion
#'   printString realpath rowData rowData<- rowRanges rowRanges<- sampleData
#'   sampleNames separator showSlotInfo standardizeCall
#' @importFrom bcbioBase getBarcodeCutoffFromCommands getGTFFileFromYAML
#'   getLevelFromCommands getSampleDataFromYAML getUMITypeFromCommands
#'   importDataVersions importProgramVersions projectDir runDate sampleDirs
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
#' @importFrom stringr str_extract
"_PACKAGE"



## REWORK AS ACIDCLI.
#' @importFrom cli cat_line cli_alert cli_alert_success cli_alert_warning cli_h1
#'   cli_h2 cli_text
NULL
