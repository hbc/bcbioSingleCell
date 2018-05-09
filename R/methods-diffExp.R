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
#' @param zeroWeights Package to use for zero weight calculations. Defaults to
#'   zinbwave but zingeR is also supported.
#' @param caller Package to use for differential expression calling. Defaults
#'   to edgeR (faster for large datasets) but DESeq2 is also supported.
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
#' data <- metrics(cellranger_small)
#' numerator <- rownames(data)[which(data[["sampleName"]] == "proximal")]
#' denominator <- rownames(data)[which(data[["sampleName"]] == "distal")]
#'
#' # zingeR-edgeR
#' x <- diffExp(
#'     object = cellranger_small,
#'     numerator = numerator,
#'     denominator = denominator,
#'     zeroWeights = "zingeR",
#'     caller = "edgeR"
#' )
#' class(x)
#' glimpse(x[["table"]])
#'
#' # seurat ====
#' # Expression in cluster 3 relative to cluster 2
#' numerator <- Seurat::WhichCells(Seurat::pbmc_small, ident = 3L)
#' denominator <- Seurat::WhichCells(Seurat::pbmc_small, ident = 2L)
#'
#' # zingeR-edgeR
#' x <- diffExp(
#'     object = Seurat::pbmc_small,
#'     numerator = numerator,
#'     denominator = denominator,
#'     zeroWeights = "zingeR",
#'     caller = "edgeR"
#' )
#' class(x)
#' glimpse(x[["table"]])
NULL



maxit <- 1000L



# zingeR + edgeR analysis
# Note that TMM needs to be consistently applied for both
# `calcNormFactors()` and `zeroWeightsLS()`
.diffExp.zingeR.edgeR <- function(  # nolint
    counts,
    design,
    group
) {
    message("zingeR-edgeR")
    counts <- as.matrix(counts)
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



.diffExp.zingeR.DESeq2 <- function(  # nolint
    counts,
    design,
    group
) {
    message("zingeR-DESeq2")
    counts <- as.matrix(counts)
    colData <- data.frame(group = group)
    designFormula <- formula("~ group")
    dds <- DESeq2::DESeqDataSetFromMatrix(
        countData = counts,
        colData = colData,
        design = designFormula
    )
    weights <- zingeR::zeroWeightsLS(
        counts = counts(dds),
        design = design,
        maxit = maxit,
        normalization = "DESeq2_poscounts",
        colData = colData,
        designFormula = designFormula
    )
    assays(dds)[["weights"]] <- weights
    dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq2::estimateDispersions(dds)
    dds <- DESeq2::nbinomWaldTest(
        dds,
        betaPrior = TRUE,
        useT = TRUE,
        df = rowSums(weights) - 2L
    )
    res <- DESeq2::results(dds)
    res
}



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
        zeroWeights = c("zingeR", "zinbwave"),
        caller = c("edgeR", "DESeq2")
    ) {
        assert_is_character(numerator)
        assert_is_character(denominator)
        assert_are_disjoint_sets(numerator, denominator)
        zeroWeights <- match.arg(zeroWeights)
        requireNamespace(zeroWeights)
        caller <- match.arg(caller)
        requireNamespace(caller)

        # Counts matrix
        cells <- c(numerator, denominator)
        counts <- assay(object)
        counts <- as.matrix(counts)
        counts <- counts[, cells]
        assert_has_dimnames(counts)

        message(paste(
            "Low gene count filter",
            "Requiring at least 25 cells with counts of 5 or more",
            sep = "\n"
        ))
        keep <- rowSums(counts >= 5L) >= 25L
        print(table(keep))
        counts <- counts[keep, ]

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

        fun <- get(paste(
            "",
            "diffExp",
            zeroWeights,
            caller,
            sep = "."
        ))
        fun(
            counts = counts,
            design = design,
            group = group
        )
    }
)



#' @rdname diffExp
#' @export
setMethod(
    "diffExp",
    signature("seurat"),
    getMethod("diffExp", "SingleCellExperiment")
)
