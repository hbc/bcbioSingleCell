#' Fetch Dimensionality Reduction Data
#'
#' Used by the [fetchPCAData()] and [fetchTSNEData()] functions.
#'
#' @author Rory Kirchner, Michael Steinbaugh
#' @keywords internal
#' @noRd
#'
#' @importFrom basejump camel
#' @importFrom dplyr group_by mutate mutate_if ungroup
#' @importFrom Seurat FetchData
#' @importFrom tibble rownames_to_column
#' @importFrom rlang !! sym
#'
#' @param object `seurat`.
#' @param dimCode Character vector of X and Y coordinate data to be used for
#'   plotting. This can be `c("tSNE_1", "tSNE_2")` for tSNE data, or `c("PC1",
#'   "PC2")` for PCA data.
#'
#' @return `data.frame`.
#'
#' @examples
#' .fetchDRData.seurat(
#'     pbmc_small,
#'     dimCode = c(x = "tSNE_1", y = "tSNE_2")
#'  ) %>%
#'     glimpse()
.fetchDRData.seurat <- function(object, dimCode) {  # nolint
    fetch <- Seurat::FetchData(object, vars.all = dimCode)
    ident <- slot(object, "ident")
    metadata <- slot(object, "meta.data")
    assert_are_identical(rownames(fetch), names(ident))
    assert_are_identical(rownames(fetch), rownames(metadata))

    # Bind into a single data.frame
    data <- cbind(
        fetch,
        as.data.frame(ident),
        metadata
    ) %>%
        # Note that this step may make the resolution columns confusing
        # (e.g. `res08` instead of `res.0.8`)
        camel() %>%
        # Move rownames before performing tidyverse operations
        rownames_to_column() %>%
        .sanitizeMetrics() %>%
        # Group by ident here for center calculations
        group_by(!!sym("ident")) %>%
        mutate(
            centerX = median(.data[[camel(dimCode[[1L]], strict = FALSE)]]),
            centerY = median(.data[[camel(dimCode[[2L]], strict = FALSE)]])
        ) %>%
        ungroup() %>%
        as.data.frame() %>%
        column_to_rownames()
    data
}
