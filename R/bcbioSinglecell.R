#' bcbioSinglecell
#'
#' Utility functions for analysis of bcbio-nextgen single-cell RNA-seq data.
#'
#' Consult the
#' \href{http://bioinformatics.sph.harvard.edu/bcbioSinglecell}{package website}
#' for more information.
#'
#' @docType package
#' @name bcbioSinglecell
#' @keywords internal
NULL



# Globals ====
globalVariables(".")
bins <- 100
fail_color <- "red"
pass_color <- "green"
warn_color <- "orange"



# Imports ====
## General ----
#' @importFrom biomaRt getBM listMarts useEnsembl
#' @importFrom knitr asis_output kable opts_knit
#' @importFrom methods as show
#' @importFrom R.utils gunzip gzip
#' @importFrom stats aggregate median
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils globalVariables sessionInfo
NULL

## Single-cell ----
## @import scater
## Namespace collisons with dplyr: arrange, filter, mutate, rename

## Sparse matrices ----
#' @importFrom Matrix colSums readMM writeMM
#' @importFrom Matrix.utils aggregate.Matrix
NULL

## Visualization ----
#' @importFrom graphics hist
#' @importFrom scales math_format trans_breaks trans_format
NULL

## tidyverse ----
## http://tidyverse.org/
#' @import dplyr
#' @import ggplot2
#' @import readr
#' @import readxl
#' @import stringr
#' @importFrom magrittr %>% set_names set_rownames
#' @importFrom rlang !! .data sym
#' @importFrom tibble as_tibble glimpse remove_rownames tibble
#' @importFrom tidyr expand_ separate_
NULL



# Re-exports ====
#' @usage NULL
#' @export
magrittr::`%>%`

#' @usage NULL
#' @export
dplyr::arrange

#' @usage NULL
#' @export
Matrix::colSums

#' @usage NULL
#' @export
dplyr::filter

#' @usage NULL
#' @export
tibble::glimpse

#' @usage NULL
#' @export
dplyr::group_by

#' @usage NULL
#' @export
dplyr::left_join

#' @usage NULL
#' @export
dplyr::mutate

#' @usage NULL
#' @export
readxl::read_excel

#' @usage NULL
#' @export
tibble::remove_rownames

#' @usage NULL
#' @export
dplyr::select

#' @usage NULL
#' @export
dplyr::top_n

#' @usage NULL
#' @export
readr::write_csv
