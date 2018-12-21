#' bcbioSingleCell
#'
#' Import and analyze [bcbio](https://bcbio-nextgen.readthedocs.io/) single-cell
#' RNA-seq data.
#'
#' @aliases NULL
#' @keywords internal
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @importMethodsFrom basejump coerce
#'
#' @importFrom Biobase sampleNames
#' @importFrom BiocGenerics cbind counts do.call
#' @importFrom DropletUtils barcodeRanks
#' @importFrom Matrix colSums readMM rowSums
#' @importFrom S4Vectors as.data.frame as.matrix mcols mcols<- merge metadata
#'   metadata<-
#' @importFrom SingleCellExperiment isSpike
#' @importFrom SummarizedExperiment assayNames assays colData colData<-
#'   rowRanges rowRanges<-
#' @importFrom basejump basejump_geom_abline basejump_geom_label
#'   basejump_geom_label_average basejump_geom_label_repel camel cell2sample
#'   detectLanes emptyRanges formalsList import interestingGroups
#'   interestingGroups<- makeDimnames makeGRangesFromEnsembl makeGRangesFromGFF
#'   makeNames makeSingleCellExperiment mapCellsToSamples markdownHeader
#'   markdownPlots matchArgsToDoCall matchInterestingGroups metrics
#'   minimalSampleData plotZerosVsDepth prepareTemplate printString realpath
#'   sampleData separator showSlotInfo
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
#' @importFrom goalie areDisjointSets areSetEqual assert hasLength hasNames
#'   hasRows hasRownames isADirectory isAFile isAURL isAny isCharacter
#'   isDirectory isFile isFlag isGGScale isHeaderLevel isInLeftOpenRange
#'   isInRange isInRightOpenRange isInt isNonEmpty isNonNegative isNumber
#'   isPositive isString isSubset validate validateClasses
#' @importFrom graphics hist
#' @importFrom magrittr %>%
#' @importFrom methods .hasSlot as as<- is new show slot slot<- validObject
#' @importFrom readr read_lines read_tsv
#' @importFrom rlang !! := sym syms UQ
#' @importFrom scales percent
#' @importFrom stats ecdf
#' @importFrom stringr str_extract str_pad
#' @importFrom tibble as_tibble tibble
#' @importFrom utils capture.output globalVariables packageVersion
"_PACKAGE"
