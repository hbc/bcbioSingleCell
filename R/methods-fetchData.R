#' Fetch Data Functions
#'
#' @name fetchData
#' @family Data Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return
#' - `fetchGeneData()`: `matrix`.
#' - `fetchTSNEExpressionData()`: `data.frame` containing t-SNE coordinates,
#'   sample metadata, and aggregate marker expression values (`mean`, `median`,
#'   and `sum`).
#' - Other functions: `data.frame`, containing metrics.
#'
#' @examples
#' # SingleCellExperiment ====
#' object <- indrops_small
#' genes <- head(rownames(object))
#'
#' # fetchGeneData
#' x <- fetchGeneData(object, genes = genes)
#' glimpse(x)
#'
#' # seurat ====
#' object <- seurat_small
#' genes <- head(rownames(object))
#'
#' # t-SNE
#' x <- fetchTSNEData(object)
#' glimpse(x)
#'
#' # PCA
#' x <- fetchPCAData(object)
#' glimpse(x)
#'
#' # UMAP
#' x <- fetchUMAPData(object)
#' glimpse(x)
#'
#' # t-SNE gene expression
#' x <- fetchTSNEExpressionData(seurat_small, genes = genes)
#' glimpse(x)
#'
#' # UMAP gene expession
#' genes <- head(rownames(seurat_small))
#' x <- fetchUMAPExpressionData(seurat_small, genes = genes)
#' glimpse(x)
NULL



# Constructors =================================================================
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
        column_to_rownames()
}



# Methods ======================================================================
#' @rdname fetchData
#' @export
setMethod(
    "fetchGeneData",
    signature("SingleCellExperiment"),
    function(object, genes) {
        assert_is_subset(genes, rownames(object))
        counts(object) %>%
            .[genes, , drop = FALSE] %>%
            as.matrix() %>%
            t()
    }
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchGeneData",
    signature("seurat"),
    getMethod("fetchGeneData", "SingleCellExperiment")
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchPCAData",
    signature("seurat"),
    function(object) {
        .fetchDimensionalReduction.seurat(object, "pca")
    }
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchTSNEData",
    signature("seurat"),
    function(object) {
        .fetchDimensionalReduction.seurat(object, "tsne")
    }
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchTSNEExpressionData",
    signature("seurat"),
    function(object, genes) {
        assert_is_subset(genes, rownames(object))
        tsne <- fetchTSNEData(object)
        data <- fetchGeneData(object, genes = genes)
        mean <- rowMeans(data)
        median <- rowMedians(data)
        sum <- rowSums(data)
        cbind(tsne, mean, median, sum)
    }
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchUMAPData",
    signature("seurat"),
    function(object) {
        .fetchDimensionalReduction.seurat(object, "umap")
    }
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchUMAPExpressionData",
    signature("seurat"),
    function(object, genes) {
        assert_is_subset(genes, rownames(object))
        umap <- fetchUMAPData(object)
        data <- fetchGeneData(object, genes = genes)
        mean <- rowMeans(data)
        median <- rowMedians(data)
        sum <- rowSums(data)
        cbind(umap, mean, median, sum)
    }
)
