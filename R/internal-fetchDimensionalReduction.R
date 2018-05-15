#' Fetch Dimensionality Reduction Data
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
#' # t-SNE
#' x <- .fetchDimensionalReduction.seurat(
#'     object = seurat_small,
#'     reductionType = "tsne"
#' )
#' glimpse(x)
#'
#' # PCA
#' x <- .fetchDimensionalReduction.seurat(
#'     object = seurat_small,
#'     reductionType = "pca"
#' )
#' glimpse(x)
#'
#' # UMAP
#' x <- .fetchDimensionalReduction.seurat(
#'     object = seurat_small,
#'     reductionType = "umap"
#' )
#' glimpse(x)
.fetchDimensionalReduction.seurat <- function(  # nolint
    object,
    reductionType = c("tsne", "umap", "pca")
) {
    reductionType <- match.arg(reductionType)
    data <- GetDimReduction(
        object = object,
        reduction.type = reductionType,
        slot = "cell.embeddings"
    )
    # Limit to the first two columns.
    # PCA returns multiple columns, for example.
    data <- data[, seq_len(2L)]
    dimCols <- colnames(data)
    metrics <- metrics(object)
    assert_are_identical(rownames(data), rownames(metrics))
    cbind(metrics, data) %>%
        rownames_to_column() %>%
        # Group by ident here for center calculations
        group_by(!!sym("ident")) %>%
        mutate(
            centerX = median(!!sym(dimCols[[1L]])),
            centerY = median(!!sym(dimCols[[2L]]))
        ) %>%
        ungroup() %>%
        as.data.frame() %>%
        column_to_rownames() %>%
        camel()
}
