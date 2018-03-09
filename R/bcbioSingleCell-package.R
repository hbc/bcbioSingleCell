#' bcbioSingleCell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @rdname bcbioSingleCell-package
#' @name bcbioSingleCell-package
#' @docType package
#'
#' @import methods
#' @importClassesFrom Seurat seurat
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom assertive assert_all_are_dirs
#' @importFrom assertive assert_all_are_existing_files
#' @importFrom assertive assert_all_are_greater_than_or_equal_to
#' @importFrom assertive assert_all_are_matching_regex
#' @importFrom assertive assert_all_are_positive
#' @importFrom assertive assert_any_are_matching_regex
#' @importFrom assertive assert_are_identical
#' @importFrom assertive assert_has_no_duplicates
#' @importFrom assertive assert_is_a_bool
#' @importFrom assertive assert_is_a_string
#' @importFrom assertive assert_is_an_integer
#' @importFrom assertive assert_is_any_of
#' @importFrom assertive assert_is_character
#' @importFrom assertive assert_is_data.frame
#' @importFrom assertive assert_is_environment
#' @importFrom assertive assert_is_factor
#' @importFrom assertive assert_is_non_empty
#' @importFrom assertive assert_is_of_length
#' @importFrom assertive assert_is_subset
#' @importFrom assertive has_dims
#' @importFrom assertive has_names
#' @importFrom assertive is_a_string
#' @importFrom assertive is_character
#' @importFrom basejump assertIsAStringOrNULL
#' @importFrom basejump assertIsAnImplicitInteger
#' @importFrom basejump assertIsAnImplicitIntegerOrNULL
#' @importFrom basejump assertIsGene2symbol
#' @importFrom bcbioBase assertFormalInterestingGroups
#' @importFrom magrittr %>%
#' @importFrom rlang !! .data abort inform sym warn
#' @importFrom S4Vectors aggregate metadata
#' @importFrom SummarizedExperiment assay assays colData rowData
NULL

#' @importFrom utils globalVariables
globalVariables(".")

#' @importFrom utils packageVersion
packageVersion <- packageVersion("bcbioSingleCell")

# Trailing number is to match cellranger output
barcodePattern <- ")_([ACGT_]{6,})(_[0-9]+)?$"
lanePattern <- "_L(\\d{3})"
metadataPriorityCols <- c("sampleID", "sampleName", "description")
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"

sepBar <- "============================================================"
