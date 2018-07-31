#' Fetch Data Functions
#'
#' @note Some of these functions require "`ident`" to be defined in [colData()].
#'   For t-SNE, UMAP, and PCA, the reduced dimensions must be defined in the
#'   [reducedDims()] slot.
#'
#' @name fetchData
#' @family Data Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#' @param minimal `boolean`. Return minimal data without metrics.
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
#' # Genes
#' x <- fetchGeneData(object, genes = genes)
#' glimpse(x)
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
#' x <- fetchTSNEExpressionData(object, genes = genes)
#' glimpse(x)
#'
#' # UMAP gene expession
#' x <- fetchUMAPExpressionData(object, genes = genes)
#' glimpse(x)
NULL



# Constructors =================================================================
.assertHasIdent <- function(object) {
    stopifnot(is(object, "SingleCellExperiment"))
    assert_is_subset("ident", colnames(colData(object)))
}



.reducedDimsData <- function(
    object,
    reduction = c("PCA", "TSNE", "UMAP"),
    minimal = FALSE
) {
    object <- as(object, "SingleCellExperiment")
    .assertHasIdent(object)
    reduction <- match.arg(reduction)
    assert_is_a_bool(minimal)

    data <- slot(object, "reducedDims")[[reduction]]
    if (!is.matrix(data)) {
        stop(
            paste(reduction, "dimensional reduction has not been calculated"),
            call. = FALSE
        )
    }
    if (isTRUE(minimal)) {
        return(data)
    }

    # Limit to the first two columns. PCA returns multiple columns.
    data <- as.data.frame(data)[, seq_len(2L)]
    dimCols <- colnames(data)
    colData <- as.data.frame(colData(object))
    assert_are_identical(rownames(data), rownames(colData))
    cbind(colData, data) %>%
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
    signature("SingleCellExperiment"),
    function(object, minimal = FALSE) {
        .reducedDimsData(
            object = object,
            reduction = "PCA",
            minimal = minimal
        )
    }
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchPCAData",
    signature("seurat"),
    getMethod("fetchPCAData", "SingleCellExperiment")
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchTSNEData",
    signature("SingleCellExperiment"),
    function(object, minimal = FALSE) {
        .reducedDimsData(
            object = object,
            reduction = "TSNE",
            minimal = minimal
        )
    }
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchTSNEData",
    signature("seurat"),
    getMethod("fetchTSNEData", "SingleCellExperiment")
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchTSNEExpressionData",
    signature("SingleCellExperiment"),
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
    "fetchTSNEExpressionData",
    signature("seurat"),
    getMethod("fetchTSNEExpressionData", "SingleCellExperiment")
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchUMAPData",
    signature("SingleCellExperiment"),
    function(object, minimal = FALSE) {
        .reducedDimsData(
            object,
            reduction = "UMAP",
            minimal = minimal
        )
    }
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchUMAPData",
    signature("seurat"),
    getMethod("fetchUMAPData", "SingleCellExperiment")
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchUMAPExpressionData",
    signature("SingleCellExperiment"),
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



#' @rdname fetchData
#' @export
setMethod(
    "fetchUMAPExpressionData",
    signature("seurat"),
    getMethod("fetchUMAPExpressionData", "SingleCellExperiment")
)
