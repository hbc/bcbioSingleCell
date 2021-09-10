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
#'   acid_geom_label_average acid_geom_label_repel autoDiscreteColorScale
#'   autoDiscreteFillScale sym syms
#' @importFrom BiocParallel bplapply bpmapply bpparam
#' @importFrom basejump DataFrame DataFrameList SimpleList abort alert
#'   alertSuccess alertWarning assayNames assay assays assays<- as_tibble
#'   calculateMetrics camelCase capture.output cbind cell2sample colData
#'   colData<- counts detectLanes do.call droplevels emptyRanges formalsList h1
#'   h2 import importSampleData interestingGroups interestingGroups<- lapply
#'   leftJoin makeDimnames makeGRangesFromEnsembl makeGRangesFromGFF makeLabel
#'   makeNames makeSingleCellExperiment mapCellsToSamples markdownPlots
#'   matchInterestingGroups mcols mcols<- metadata metadata<- metrics
#'   metricsCols minimalSampleData packageName packageVersion printString
#'   realpath rowData rowData<- rowRanges rowRanges<- sampleData sampleNames
#'   separator showSlotInfo standardizeCall str_extract toInlineString
#' @importFrom bcbioBase getBarcodeCutoffFromCommands getGTFFileFromYAML
#'   getLevelFromCommands getSampleDataFromYAML getUMITypeFromCommands
#'   importDataVersions importProgramVersions projectDir runDate sampleDirs
#' @importFrom ggplot2 aes facet_wrap geom_boxplot geom_histogram geom_step
#'   geom_violin ggplot labs scale_x_continuous scale_y_continuous stat_ecdf
#' @importFrom ggridges geom_density_ridges
#' @importFrom goalie allAreDirectories allAreFiles areDisjointSets areSetEqual
#'   assert hasLength hasNames hasRownames hasValidDimnames isADirectory isAFile
#'   isAURL isAny isCharacter isDirectory isFile isFlag isInt isString isSubset
#'   validate validateClasses
#' @importFrom graphics hist
#' @importFrom methods as as<- is new setClass show slot slot<- validObject
#'   .hasSlot
"_PACKAGE"
