#' bcbioSingleCell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
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
#' @importFrom assertive.base assert_are_identical
#' @importFrom assertive.files assert_all_are_dirs assert_all_are_existing_files
#' @importFrom assertive.numbers assert_all_are_in_left_open_range
#'   assert_all_are_in_range assert_all_are_in_right_open_range
#'   assert_all_are_non_negative assert_all_are_positive
#' @importFrom assertive.properties assert_has_names assert_has_rows
#'   assert_is_non_empty
#' @importFrom assertive.sets assert_are_disjoint_sets assert_are_set_equal
#'   assert_is_subset
#' @importFrom assertive.types assert_is_a_bool assert_is_a_number
#'   assert_is_a_string assert_is_all_of assert_is_an_integer assert_is_any_of
#'   assert_is_character assert_is_data.frame assert_is_factor assert_is_integer
#'   assert_is_list assert_is_numeric is_a_string
#' @importFrom assertthat assert_that validate_that
#' @importFrom basejump emptyRanges makeGRangesFromEnsembl makeGRangesFromGFF
#' @importFrom basejump basejump_geom_abline basejump_geom_label
#'   basejump_geom_label_average basejump_geom_label_repel camel cell2sample
#'   detectLanes formalsList import interestingGroups interestingGroups<-
#'   makeDimnames makeNames makeSingleCellExperiment mapCellsToSamples
#'   markdownHeader markdownPlotlist matchArgsToDoCall matchInterestingGroups
#'   metrics minimalSampleData plotZerosVsDepth prepareTemplate printString
#'   realpath sampleData separator showSlotInfo
#' @importFrom bcbioBase getBarcodeCutoffFromCommands getLevelFromCommands
#'   getSampleDataFromYAML getUMITypeFromCommands projectDir readDataVersions
#'   readProgramVersions readSampleData runDate sampleDirs
#' @importFrom cowplot plot_grid
#' @importFrom dplyr arrange bind_rows desc filter group_by left_join mutate
#'   mutate_all mutate_if n pull slice summarise
#' @importFrom ggplot2 aes facet_wrap geom_bar geom_boxplot geom_histogram
#'   geom_hline geom_line geom_point geom_smooth geom_step geom_violin ggplot
#'   labs scale_x_continuous scale_y_continuous stat_ecdf theme
#' @importFrom ggridges geom_density_ridges
#' @importFrom goalie assertHasRownames assertIsAnImplicitInteger
#'   assertIsAnImplicitIntegerOrNULL assertIsColorScaleDiscreteOrNULL
#'   assertIsFillScaleDiscreteOrNULL assertIsHeaderLevel assertIsStringOrNULL
#'   hasRownames
#' @importFrom graphics hist
#' @importFrom magrittr %>%
#' @importFrom methods .hasSlot as as<- is new show slot slot<- validObject
#' @importFrom readr read_lines read_tsv
#' @importFrom rlang !! := has_length sym syms UQ
#' @importFrom scales percent
#' @importFrom stats ecdf
#' @importFrom stringr str_extract str_pad
#' @importFrom tibble as_tibble tibble
#' @importFrom utils capture.output globalVariables packageVersion
"_PACKAGE"
