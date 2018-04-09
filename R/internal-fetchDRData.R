#' Fetch Dimensionality Reduction Data
#'
#' Used by the [fetchPCAData()] and [fetchTSNEData()] functions.
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @inheritParams general
#' @param dimCode Character vector of X and Y coordinate data to be used for
#'   plotting. This can be `c("tSNE_1", "tSNE_2")` for tSNE data, or `c("PC1",
#'   "PC2")` for PCA data.
#'
#' @return `data.frame`.
#'
#' @examples
#' x <- .fetchDRData.seurat(
#'     object = Seurat::pbmc_small,
#'     dimCode = c(x = "tSNE_1", y = "tSNE_2")
#' )
#' glimpse(x)
.fetchDRData.seurat <- function(object, dimCode) {  # nolint
    fetch <- Seurat::FetchData(object, vars.all = dimCode)
    metrics <- metrics(object)
    assert_are_identical(rownames(fetch), rownames(metrics))
    cbind(fetch, metrics) %>%
        rownames_to_column() %>%
        camel() %>%
        # Group by ident here for center calculations
        group_by(!!sym("ident")) %>%
        mutate(
            centerX = median(.data[[camel(dimCode[[1L]], strict = FALSE)]]),
            centerY = median(.data[[camel(dimCode[[2L]], strict = FALSE)]])
        ) %>%
        ungroup() %>%
        as.data.frame() %>%
        column_to_rownames()
}
