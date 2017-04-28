#' \code{bcbioSinglecell} package
#'
#' Utility functions for analysis of \code{bcbio-nextgen} single-cell RNA-seq
#' data.
#'
#' See the README on \href{https://github.com/roryk/bcbioSinglecell}{GitHub}.
#'
#' @docType package
#' @name bcbioSinglecell
NULL



# Globals ====
globalVariables(".")



# Imports ====
## General ----
#' @importFrom biomaRt getBM useEnsembl
#' @importFrom knitr asis_output kable opts_knit
#' @importFrom methods as show
#' @importFrom R.utils gzip
#' @importFrom stats aggregate median
#' @importFrom tools file_ext file_path_sans_ext
#' @importFrom utils globalVariables sessionInfo
NULL

## Sparse matrices ----
#' @importFrom Matrix readMM writeMM
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
#' @importFrom tibble as_tibble glimpse tibble
#' @importFrom tidyr expand_ separate_
NULL



# Re-exports ====
#' @usage NULL
#' @export
magrittr::`%>%`

#' @usage NULL
#' @export
tibble::glimpse
