## S4 classes ==================================================================

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
NULL



## S4 generics and methods =====================================================

#' @importFrom AcidExperiment sampleNames
#' @importFrom AcidGenerics calculateMetrics camelCase cell2sample
#' interestingGroups interestingGroups<- leftJoin makeDimnames makeLabel
#' makeNames metrics plotReadsPerCell sampleData
#' @importFrom BiocGenerics counts updateObject
#' @importFrom BiocParallel bplapply bpmapply
#' @importFrom S4Vectors cbind do.call droplevels lapply mcols mcols<-
#' metadata metadata<-
#' @importFrom SummarizedExperiment assayNames assay assays assays<- colData
#' colData<- rowData rowData<- rowRanges rowRanges<-
#' @importFrom methods coerce show
#' @importFrom pipette import
#'
#' @importMethodsFrom AcidExperiment calculateMetrics interestingGroups
#' interestingGroups<- metrics sampleData sampleNames
#' @importMethodsFrom AcidPlyr leftJoin
#' @importMethodsFrom AcidSingleCell cell2sample sampleData
#' @importMethodsFrom pipette coerce import
#' @importMethodsFrom syntactic camelCase makeDimnames makeLabel makeNames
NULL



## S3 generics =================================================================

#' @importFrom pipette as_tibble
NULL



## Standard functions ==========================================================

#' @importFrom AcidBase metricsCols printString realpath showSlotInfo
#' standardizeCall
#' @importFrom AcidCLI abort alert alertSuccess alertWarning h1 h2
#' separator toInlineString
#' @importFrom AcidExperiment detectLanes importSampleData
#' matchInterestingGroups minimalSampleData
#' @importFrom AcidGenomes emptyRanges makeGRangesFromEnsembl makeGRangesFromGFF
#' @importFrom AcidMarkdown markdownPlots
#' @importFrom AcidPlots !! acid_geom_abline acid_geom_label
#' acid_geom_label_average acid_geom_label_repel autoDiscreteColorScale
#' autoDiscreteFillScale sym syms
#' @importFrom AcidSingleCell makeSingleCellExperiment mapCellsToSamples
#' @importFrom IRanges DataFrameList
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom BiocParallel bpparam
#' @importFrom bcbioBase getBarcodeCutoffFromCommands getGTFFileFromYAML
#' getLevelFromCommands getSampleDataFromYAML getUMITypeFromCommands
#' importDataVersions importProgramVersions projectDir runDate sampleDirs
#' @importFrom ggplot2 aes facet_wrap geom_boxplot geom_histogram geom_step
#' geom_violin ggplot labs scale_x_continuous scale_y_continuous stat_ecdf
#' @importFrom ggridges geom_density_ridges
#' @importFrom goalie allAreDirectories allAreFiles areDisjointSets areSetEqual
#' assert hasLength hasNames hasRownames hasValidDimnames isADirectory isAFile
#' isAURL isAny isCharacter isDirectory isFile isFlag isInt isString isSubset
#' validate validateClasses
#' @importFrom graphics hist
#' @importFrom methods as as<- is new setClass slot slot<- validObject
#' .hasSlot
#' @importFrom stringr str_extract
#' @importFrom utils capture.output packageName packageVersion
NULL
