#' bcbioSingleCell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @name bcbioSingleCell-package
#' @docType package
#'
#' @importClassesFrom Seurat seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Biobase rowMedians
#' @importFrom GenomicFeatures genes makeTxDbFromGFF transcripts
#' @importFrom Matrix colSums readMM rowMeans rowSums
#' @importFrom Matrix.utils aggregate.Matrix
#' @importFrom S4Vectors DataFrame aggregate as.data.frame as.matrix
#'   complete.cases mcols metadata metadata<- na.omit
#' @importFrom SingleCellExperiment SingleCellExperiment isSpike isSpike<-
#'   spikeNames
#' @importFrom SummarizedExperiment assay assays colData rowData rowRanges
#' @importFrom basejump assignAndSaveData camel convertGenesToSymbols
#'   convertUCSCBuildToEnsembl detectOrganism dynamicPlotlist emptyRanges
#'   initializeDirectory makeGRangesFromEnsembl makeGRangesFromGFF makeNames
#'   makeTx2geneFromGFF markdownHeader readFileByExtension readYAML
#'   sanitizeSampleData theme_midnight theme_paperwhite tx2geneFromGFF
#' @importFrom bcbioBase flatFiles prepareSummarizedExperiment readDataVersions
#'   readLog readProgramVersions readSampleData readTx2gene
#'   sampleYAMLMetadata uniteInterestingGroups
#' @importFrom cowplot draw_plot ggdraw plot_grid
#' @importFrom dplyr arrange bind_rows desc filter group_by group_vars left_join
#'   matches mutate mutate_all mutate_if n select select_if summarize
#'   summarize_all ungroup
#' @importFrom ggplot2 aes_ aes_string coord_flip element_blank element_line
#'   element_rect element_text expand_limits facet_wrap geom_bar geom_boxplot
#'   geom_histogram geom_hline geom_label geom_line geom_point geom_smooth
#'   geom_step geom_text geom_violin geom_vline ggtitle guide_colorbar
#'   guide_legend guides labs qplot scale_radius scale_x_continuous
#'   scale_y_continuous stat_ecdf theme xlab xlim ylab
#' @importFrom ggridges geom_density_ridges
#' @importFrom graphics hist
#' @importFrom grid unit
#' @importFrom jsonlite read_json
#' @importFrom knitr kable
#' @importFrom magrittr %>% set_colnames set_names set_rownames
#' @importFrom methods .hasSlot as as<- getMethod is new show slot slot<-
#'   validObject
#' @importFrom parallel mclapply mcmapply
#' @importFrom pbapply pblapply
#' @importFrom readr read_lines read_tsv
#' @importFrom rlang !! !!! sym syms UQ
#' @importFrom scales percent pretty_breaks
#' @importFrom stats as.formula ecdf fitted median model.matrix predict relevel
#'   reorder smooth.spline
#' @importFrom stringr str_extract str_match str_pad str_split
#' @importFrom tibble as_tibble column_to_rownames has_rownames remove_rownames
#'   rownames_to_column tibble
#' @importFrom tidyr gather
#' @importFrom utils globalVariables packageVersion
#'
#' @importFrom assertive.base assert_are_identical
#' @importFrom assertive.files assert_all_are_dirs
#' @importFrom assertive.files assert_all_are_existing_files
#' @importFrom assertive.numbers assert_all_are_greater_than_or_equal_to
#' @importFrom assertive.numbers assert_all_are_in_left_open_range
#' @importFrom assertive.numbers assert_all_are_in_range
#' @importFrom assertive.numbers assert_all_are_in_right_open_range
#' @importFrom assertive.numbers assert_all_are_non_negative
#' @importFrom assertive.numbers assert_all_are_positive
#' @importFrom assertive.properties assert_are_same_length
#' @importFrom assertive.properties assert_has_dimnames
#' @importFrom assertive.properties assert_has_names
#' @importFrom assertive.properties assert_has_no_duplicates
#' @importFrom assertive.properties assert_has_rows
#' @importFrom assertive.properties assert_is_non_empty
#' @importFrom assertive.properties assert_is_of_length
#' @importFrom assertive.properties has_dims
#' @importFrom assertive.properties has_names
#' @importFrom assertive.sets assert_are_disjoint_sets
#' @importFrom assertive.sets assert_are_intersecting_sets
#' @importFrom assertive.sets assert_is_subset
#' @importFrom assertive.strings assert_all_are_matching_regex
#' @importFrom assertive.strings assert_all_are_non_missing_nor_empty_character
#' @importFrom assertive.strings assert_any_are_matching_regex
#' @importFrom assertive.types assert_is_a_bool
#' @importFrom assertive.types assert_is_a_number
#' @importFrom assertive.types assert_is_a_string
#' @importFrom assertive.types assert_is_all_of
#' @importFrom assertive.types assert_is_an_integer
#' @importFrom assertive.types assert_is_any_of
#' @importFrom assertive.types assert_is_character
#' @importFrom assertive.types assert_is_data.frame
#' @importFrom assertive.types assert_is_environment
#' @importFrom assertive.types assert_is_factor
#' @importFrom assertive.types assert_is_function
#' @importFrom assertive.types assert_is_integer
#' @importFrom assertive.types assert_is_list
#' @importFrom assertive.types assert_is_numeric
#' @importFrom assertive.types assert_is_tbl_df
#' @importFrom assertive.types is_a_string
#' @importFrom assertive.types is_character
#'
#' @importFrom basejump assertHasRownames
#' @importFrom basejump assertIsAHeaderLevel
#' @importFrom basejump assertIsAStringOrNULL
#' @importFrom basejump assertIsAnImplicitInteger
#' @importFrom basejump assertIsAnImplicitIntegerOrNULL
#' @importFrom basejump assertIsANumberOrNULL
#' @importFrom basejump assertIsCharacterOrNULL
#' @importFrom basejump assertIsColorScaleContinuousOrNULL
#' @importFrom basejump assertIsColorScaleDiscreteOrNULL
#' @importFrom basejump assertIsFillScaleDiscreteOrNULL
#' @importFrom basejump assertIsGene2symbol
#' @importFrom basejump assertIsTx2gene
#' @importFrom basejump hasRownames
#'
#' @importFrom bcbioBase assertFormalInterestingGroups
NULL
