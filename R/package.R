#' bcbioSingleCell
#'
#' Import and analyze [bcbio](https://bcbio-nextgen.readthedocs.io/) single-cell
#' RNA-seq data.
#'
#' @aliases NULL
#' @keywords internal
"_PACKAGE"



## S4 classes ==================================================================

#' @importClassesFrom SingleCellExperiment SingleCellExperiment
NULL



## S4 generics and methods =====================================================

#' @importFrom AcidExperiment sampleNames
#' @importFrom AcidGenerics calculateMetrics camelCase droplevels2
#' interestingGroups interestingGroups<- leftJoin makeDimnames makeLabel
#' makeNames metrics plotReadsPerCell sampleData
#' @importFrom BiocGenerics counts updateObject
#' @importFrom S4Vectors cbind do.call droplevels lapply mcols mcols<-
#' metadata metadata<-
#' @importFrom SummarizedExperiment assayNames assay assays assays<- colData
#' colData<- rowData rowData<- rowRanges rowRanges<-
#' @importFrom methods coerce show
#' @importFrom pipette import
NULL

#' @importMethodsFrom AcidExperiment calculateMetrics interestingGroups
#' interestingGroups<- metrics sampleData sampleNames
#' @importMethodsFrom AcidPlyr leftJoin
#' @importMethodsFrom AcidSingleCell sampleData
#' @importMethodsFrom pipette import
#' @importMethodsFrom syntactic camelCase makeDimnames makeLabel makeNames
NULL



## Standard functions ==========================================================

#' @importFrom AcidBase metricsCols printString realpath showSlotInfo
#' standardizeCall strMatch
#' @importFrom AcidCLI abort alert alertSuccess alertWarning h1 h2
#' separator toInlineString
#' @importFrom AcidExperiment detectLanes droplevels2 importSampleData
#' matchInterestingGroups minimalSampleData
#' @importFrom AcidGenomes emptyRanges makeGRangesFromEnsembl makeGRangesFromGFF
#' @importFrom AcidMarkdown markdownPlots
#' @importFrom AcidPlots !!! .data acid_geom_abline acid_geom_label
#' acid_geom_label_average acid_geom_label_repel acid_scale_color_discrete
#' acid_scale_fill_discrete syms
#' @importFrom AcidSingleCell makeSingleCellExperiment mapCellsToSamples
#' @importFrom IRanges DataFrameList
#' @importFrom S4Vectors DataFrame SimpleList
#' @importFrom bcbioBase getBarcodeCutoffFromCommands getGTFFileFromYAML
#' getLevelFromCommands getSampleDataFromYAML getUMITypeFromCommands
#' importDataVersions importProgramVersions projectDir runDate sampleDirs
#' @importFrom ggplot2 aes facet_wrap geom_boxplot geom_histogram geom_step
#' geom_violin ggplot labs scale_x_continuous scale_y_continuous stat_ecdf vars
#' @importFrom ggridges geom_density_ridges
#' @importFrom goalie allAreDirectories allAreFiles areDisjointSets areSetEqual
#' assert hasLength hasNames hasRownames hasValidDimnames isADirectory isAFile
#' isAURL isAny isCharacter isDirectory isFile isFlag isInt isString isSubset
#' requireNamespaces validate validateClasses
#' @importFrom methods .hasSlot as as<- is new setClass slot slot<- validObject
#' @importFrom parallel mcMap mclapply
#' @importFrom utils capture.output packageName packageVersion
NULL
