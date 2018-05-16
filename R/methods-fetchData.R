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
#' genes <- rownames(indrops_small) %>% head(2L)
#' fetchGeneData(indrops_small, genes = genes) %>% glimpse()
#'
#' # seurat ====
#' # t-SNE
#' x <- fetchTSNEData(seurat_small)
#' glimpse(x)
#'
#' # PCA
#' x <- fetchPCAData(seurat_small)
#' glimpse(x)
#'
#' # UMAP
#' x <- fetchUMAPData(seurat_small)
#' glimpse(x)
#'
#' # Gene expression (marker) t-SNE
#' genes <- head(rownames(seurat_small))
#' x <- fetchTSNEExpressionData(seurat_small, genes = genes)
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
        column_to_rownames() %>%
        camel()
}



# Methods ======================================================================
#' @rdname fetchData
#' @export
setMethod(
    "fetchGeneData",
    signature("SingleCellExperiment"),
    function(object, genes) {
        counts <- counts(object)
        assert_is_character(genes)
        genes <- make.names(genes)
        assert_are_intersecting_sets(genes, rownames(counts))
        counts[genes, , drop = FALSE] %>%
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
        assert_is_character(genes)
        genes <- make.names(genes)
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



# TODO Add `fetchUMAPExpressionData()` support
