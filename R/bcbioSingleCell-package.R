#' bcbioSingleCell
#'
#' Import and analyze [bcbio](http://bcbio-nextgen.readthedocs.io) single-cell
#' RNA-seq data.
#'
#' @rdname bcbioSingleCell-package
#' @name bcbioSingleCell-package
#' @docType package
#'
#' @import methods SummarizedExperiment
#'
#' @importClassesFrom Seurat seurat
#' @importClassesFrom SummarizedExperiment SummarizedExperiment
#'
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom S4Vectors aggregate metadata
#' @importFrom utils globalVariables
NULL

globalVariables(".")

# Trailing number is to match cellranger output
barcodePattern <- ")_([ACGT_]{6,})(_[0-9]+)?$"
lanePattern <- "_L(\\d{3})"
projectDirPattern <- "^(\\d{4}-\\d{2}-\\d{2})_([^/]+)$"
metadataPriorityCols <- c("sampleID", "sampleName", "description")
