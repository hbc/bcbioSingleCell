#' @rdname bcbioSingleCell
#' @aliases NULL
#' @exportClass bcbioSingleCell
#' @usage NULL
bcbioSingleCell <- setClass(
    Class = "bcbioSingleCell",
    contains = "SingleCellExperiment"
)



setOldClass(Classes = c("grouped_df", "tbl_df", "tibble"))
