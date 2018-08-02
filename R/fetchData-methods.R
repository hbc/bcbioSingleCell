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
#' glimpse(genes)
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



# Methods ======================================================================
#' @rdname fetchData
#' @export
setMethod(
    "fetchGeneData",
    signature("SingleCellExperiment"),
    function(
        object,
        genes,
        gene2symbol = FALSE
    ) {
        assert_is_character(genes)
        assert_has_no_duplicates(genes)
        assert_is_a_bool(gene2symbol)

        counts <- counts(object)
        assert_is_subset(genes, rownames(object))
        counts <- counts[genes, , drop = FALSE]

        # Convert gene IDs to gene names (symbols)
        if (isTRUE(gene2symbol) && !isTRUE(.useGene2symbol(object))) {
            g2s <- gene2symbol(object)
            assertIsGene2symbol(g2s)
            g2s <- g2s[genes, , drop = FALSE]
            assert_are_identical(rownames(counts), g2s[["geneID"]])
            rownames(counts) <- make.unique(g2s[["geneName"]])
        }

        t(as.matrix(counts))
    }
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchReducedDimData",
    signature("SingleCellExperiment"),
    function(
        object,
        reducedDim,
        dimsUse = c(1L, 2L),
        interestingGroups
    ) {
        object <- as(object, "SingleCellExperiment")
        .assertHasIdent(object)
        assert_is_a_string(reducedDim)
        assertIsImplicitInteger(dimsUse)
        assert_is_of_length(dimsUse, 2L)
        interestingGroups <- .prepareInterestingGroups(
            object = object,
            interestingGroups = interestingGroups
        )

        data <- slot(object, "reducedDims")[[reducedDim]]
        if (!is.matrix(data)) {
            stop(
                paste(reducedDim, "reduced dimension not calculated"),
                call. = FALSE
            )
        }
        data <- as.data.frame(data)

        metrics <- metrics(object, interestingGroups = interestingGroups)
        assert_are_identical(rownames(data), rownames(metrics))

        dimCols <- colnames(data)[dimsUse]
        assert_is_character(dimCols)

        cbind(data, metrics) %>%
            rownames_to_column() %>%
            # Group by ident here for center calculations
            group_by(!!sym("ident")) %>%
            mutate(
                x = !!sym(dimCols[[1L]]),
                y = !!sym(dimCols[[2L]]),
                centerX = median(!!sym(dimCols[[1L]])),
                centerY = median(!!sym(dimCols[[2L]]))
            ) %>%
            ungroup() %>%
            as.data.frame() %>%
            column_to_rownames()
    }
)



#' @rdname fetchData
#' @export
setMethod(
    "fetchReducedDimExpressionData",
    signature("SingleCellExperiment"),
    function(
        object,
        genes,
        reducedDim
    ) {
        assert_is_subset(genes, rownames(object))

        # Gene data
        geneData <- fetchGeneData(object = object, genes = genes)

        # Expression columns
        mean <- rowMeans(geneData)
        median <- rowMedians(geneData)
        sum <- rowSums(geneData)

        # Reduced dim data
        reducedDimData <- fetchReducedDimData(
            object = object,
            reducedDim = reducedDim
        )

        cbind(reducedDimData, mean, median, sum)
    }
)
