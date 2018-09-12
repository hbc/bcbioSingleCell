#' bcbioSingleCell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @name bcbioSingleCell-package
#' @docType package
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#'
#' @importFrom assertive.base assert_are_identical
#' @importFrom assertive.files assert_all_are_dirs assert_all_are_existing_files
#' @importFrom assertive.numbers assert_all_are_in_left_open_range
#'   assert_all_are_in_range assert_all_are_in_right_open_range
#'   assert_all_are_non_negative assert_all_are_positive
#' @importFrom assertive.properties assert_has_names assert_has_no_duplicates
#'   assert_has_rows assert_is_non_empty
#' @importFrom assertive.sets assert_are_disjoint_sets assert_are_set_equal
#'   assert_is_subset
#' @importFrom assertive.strings assert_any_are_matching_regex
#' @importFrom assertive.types assert_is_a_bool assert_is_a_number
#'   assert_is_a_string assert_is_all_of assert_is_an_integer assert_is_any_of
#'   assert_is_character assert_is_data.frame assert_is_environment
#'   assert_is_factor assert_is_function assert_is_integer assert_is_list
#'   assert_is_numeric assert_is_tbl_df is_a_string
#' @importFrom basejump assertFormalInterestingGroups assertHasRownames
#'   assertIsAHeaderLevel assertIsAStringOrNULL assertIsAnImplicitInteger
#'   assertIsAnImplicitIntegerOrNULL assertIsColorScaleDiscreteOrNULL
#'   assertIsFillScaleDiscreteOrNULL assignAndSaveData basejump_geom_abline
#'   basejump_geom_label basejump_geom_label_average basejump_geom_label_repel
#'   camel emptyRanges gene2symbol hasRownames interestingGroups
#'   interestingGroups<- initializeDirectory makeGRangesFromEnsembl
#'   makeGRangesFromGFF makeNames makeSummarizedExperiment makeTx2geneFromGFF
#'   matchCall matchInterestingGroups markdownHeader markdownPlotlist
#'   prepareTemplate printString readJSON readYAML sanitizeSampleData
#'   setArgsToDoCall stripTranscriptVersions theme_midnight theme_paperwhite
#'   tx2geneFromGFF uniteInterestingGroups
#' @importFrom bcbioBase minimalSampleData projectDir readDataVersions readLog
#'   readProgramVersions readSampleData readYAMLSampleData sampleDirs
#' @importFrom Biobase sampleNames
#' @importFrom BiocGenerics cbind counts counts<- do.call rbind
#' @importFrom BiocParallel SerialParam
#' @importFrom cowplot plot_grid
#' @importFrom crayon bold
#' @importFrom dplyr arrange bind_rows desc filter group_by left_join matches
#'   mutate mutate_all mutate_if n pull rename select select_if slice summarize
#'   summarize_all ungroup
#' @importFrom DropletUtils barcodeRanks
#' @importFrom ggplot2 aes facet_wrap geom_bar geom_boxplot geom_histogram
#'   geom_hline geom_line geom_point geom_smooth geom_step geom_text geom_violin
#'   geom_vline ggplot labs scale_x_continuous scale_y_continuous stat_ecdf
#'   theme
#' @importFrom ggridges geom_density_ridges
#' @importFrom graphics hist
#' @importFrom magrittr %>% set_colnames
#' @importFrom Matrix colSums readMM rowSums sparseMatrix
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom methods .hasSlot as as<- getMethod is new show slot slot<-
#'   validObject
#' @importFrom purrr map
#' @importFrom readr read_lines read_tsv
#' @importFrom rhdf5 h5dump h5read
#' @importFrom rlang !! !!! sym syms UQ
#' @importFrom S4Vectors DataFrame aggregate as.data.frame as.matrix
#'   complete.cases mcols mcols<- merge metadata metadata<- na.omit
#' @importFrom scales percent pretty_breaks
#' @importFrom SingleCellExperiment SingleCellExperiment isSpike isSpike<-
#'   spikeNames
#' @importFrom stats ecdf
#' @importFrom stringr str_extract str_match str_pad str_split
#' @importFrom SummarizedExperiment assay assayNames assays assays<- colData
#'   colData<- rowData rowRanges rowRanges<-
#' @importFrom tibble as_tibble column_to_rownames rownames_to_column tibble
#' @importFrom utils capture.output globalVariables packageVersion
NULL
