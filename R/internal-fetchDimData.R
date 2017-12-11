#' Fetch Data
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom basejump camel
#' @importFrom dplyr group_by mutate mutate_if ungroup
#' @importFrom Seurat FetchData
#' @importFrom stats median
#' @importFrom tibble rownames_to_column
#' @importFrom rlang !! sym
#'
#' @param object [seurat] object.
#' @param dimCode Character vector of X and Y coordinate data to be used for
#'   plotting. This can be `c("tSNE_1", "tSNE_2")` for tSNE data, or `c("PC1",
#'   "PC2")` for PCA data.
#'
#' @return [data.frame].
.fetchDimDataSeurat <- function(object, dimCode) {
    ident <- data.frame(ident = slot(object, "ident"))
    metadata <- slot(object, "meta.data")
    Seurat::FetchData(object, vars.all = dimCode) %>%
        as.data.frame() %>%
        camel(strict = FALSE) %>%
        cbind(ident, metadata) %>%
        rownames_to_column() %>%
        mutate_if(is.factor, droplevels) %>%
        # Group by ident here for center calculations
        group_by(!!sym("ident")) %>%
        mutate(
            centerX = median(.data[[camel(dimCode[[1]], strict = FALSE)]]),
            centerY = median(.data[[camel(dimCode[[2]], strict = FALSE)]])
        ) %>%
        ungroup() %>%
        as.data.frame() %>%
        column_to_rownames()
}
