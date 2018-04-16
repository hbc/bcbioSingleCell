#' Differential Expression
#'
#' @section zingeR-edgeR:
#' We first perform a differential expression analysis using zingeR posterior
#' probabilities and focussing on the count component of the ZINB model. We use
#' edgeR to model the count component.
#'
#' The weights are estimated with the core function of zingeR,
#' `zeroWeightsLS()`. It is important to be consistent with the normalization
#' procedure, i.e. if you use TMM normalization for the analysis, you should
#' also use it for estimating the zingeR posterior probabilities. The weights
#' can then be provided in the DGEList object which will be used for dispersion
#' estimation in the native edgeR `estimateDisp()` function.
#'
#' After estimation of the dispersions and zingeR posterior probabilities, the
#' `glmWeightedF()` function is used for statistical inference. This is an
#' adapted function from the `glmLRT()` function of edgeR. It uses an F-test for
#' which the denominator degrees of freedom are by default adjusted according to
#' the downweighting of excess zeros (`ZI = TRUE`). Also, independent filtering
#' can be performed on the obtained p-values (`independentFiltering = TRUE`). We
#' use the independent filtering strategy that was originally implemented in
#' DESeq2. By default, the average fitted values are used as a filter criterion.
#'
#' Consult the [zingeR vignette v2](https://goo.gl/4rTK1w) for additional
#' information.
#'
#' @name diffExp
#' @family Differential Expression Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams zingeR::zeroWeightsLS
#' @inheritParams general
#' @param numerator Group of cells to use in the numerator of the contrast
#'   (e.g. treatment).
#' @param denominator Group of cells to use in the denominator of the contrast
#'   (e.g. control).
#' @param minCells Minimum number of cells required per group.
#'
#' @seealso
#' - DESeq2: We're trying to follow the conventions used in DESeq2 for
#'   contrasts, defining the name of the factor in the design formula,
#'   numerator, and denominator level for the fold change calculations. See
#'   [DESeq2::results()] for more information.
#' - Seurat: Note that Seurat currently uses the convention `cells.1` for the
#'   numerator and `cells.2` for the denominator. See [Seurat::DiffExpTest()]
#'   for additional information.
#'
#' @return `DEGLRT`.
#'
#' @examples
#' # SingleCellExperiment ====
#' x <- metrics(cellranger_small)
#' numerator <- rownames(x)[which(x[["sampleName"]] == "proximal")]
#' denominator <- rownames(x)[which(x[["sampleName"]] == "distal")]
#' lrt <- diffExp(
#'     object = cellranger_small,
#'     numerator = numerator,
#'     denominator = denominator,
#'     maxit = 100L
#' )
#' glimpse(lrt[["table"]])
#'
#' # seurat ====
#' # Expression in cluster 3 relative to cluster 2
#' numerator <- Seurat::WhichCells(Seurat::pbmc_small, ident = 3L)
#' denominator <- Seurat::WhichCells(Seurat::pbmc_small, ident = 2L)
#' lrt <- diffExp(
#'     object = Seurat::pbmc_small,
#'     numerator = numerator,
#'     denominator = denominator,
#'     maxit = 100L
#' )
#' glimpse(lrt[["table"]])
NULL



# Methods ======================================================================
#' @rdname diffExp
#' @export
setMethod(
    "diffExp",
    signature("SingleCellExperiment"),
    function(
        object,
        numerator,
        denominator,
        minCells = 10L,
        maxit = 1000L
    ) {
        requireNamespace("zingeR")
        requireNamespace("edgeR")

        assert_is_character(numerator)
        assert_is_character(denominator)
        assert_are_disjoint_sets(numerator, denominator)
        assertIsAnImplicitInteger(minCells)
        assert_all_are_greater_than_or_equal_to(
            x = c(length(numerator), length(denominator)),
            y = minCells
        )
        assertIsAnImplicitInteger(maxit)

        # Counts matrix
        cells <- c(numerator, denominator)
        counts <- assay(object)
        counts <- counts[, cells]
        assert_has_dimnames(counts)

        # Create a cell factor to define the group for `DGEList()`
        numeratorFactor <- replicate(
            n = length(numerator),
            expr = "numerator"
        ) %>%
            factor() %>%
            set_names(numerator)
        denominatorFactor <- replicate(
            n = length(denominator),
            expr = "denominator"
        ) %>%
            factor() %>%
            set_names(denominator)
        group <- factor(c(
            as.character(numeratorFactor),
            as.character(denominatorFactor)
        ))
        names(group) <- c(names(numeratorFactor), names(denominatorFactor))
        # Ensure denominator is set as reference
        group <- relevel(group, ref = "denominator")

        # Set up the design matrix
        design <- model.matrix(~group)

        # zingeR + edgeR analysis
        # Note that TMM needs to be consistently applied for both
        # `calcNormFactors()` and `zeroWeightsLS()`
        dge <- edgeR::DGEList(counts, group = group)
        dge <- edgeR::calcNormFactors(dge, method = "TMM")

        # This is the zingeR step that is computationally expensive
        weights <- zingeR::zeroWeightsLS(
            counts = dge[["counts"]],
            design = design,
            maxit = maxit,
            normalization = "TMM"
        )

        dge[["weights"]] <- weights
        dge <- edgeR::estimateDisp(dge, design = design)
        fit <- edgeR::glmFit(dge, design = design)
        lrt <- zingeR::glmWeightedF(fit, coef = 2L, independentFiltering = TRUE)
        lrt
    }
)



#' @rdname diffExp
#' @export
setMethod(
    "diffExp",
    signature("seurat"),
    getMethod("diffExp", "SingleCellExperiment")
)
