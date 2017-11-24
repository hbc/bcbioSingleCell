#' Fetch Data
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom basejump camel
#' @importFrom dplyr left_join mutate mutate_if
#' @importFrom magrittr set_colnames set_rownames
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
    meta <- slot(object, "meta.data") %>%
        camel(strict = FALSE) %>%
        rownames_to_column("cell")
    ident <- slot(object, "ident") %>%
        as.data.frame() %>%
        set_colnames("ident") %>%
        rownames_to_column("cell")
    FetchData(object, vars.all = dimCode) %>%
        as.data.frame() %>%
        camel(strict = FALSE) %>%
        rownames_to_column("cell") %>%
        left_join(meta, by = "cell") %>%
        left_join(ident, by = "cell") %>%
        mutate_if(is.factor, droplevels) %>%
        mutate(
            centerX = median(.data[[camel(dimCode[[1]], strict = FALSE)]]),
            centerY = median(.data[[camel(dimCode[[2]], strict = FALSE)]])
        ) %>%
        set_rownames(.[["cell"]])
}
