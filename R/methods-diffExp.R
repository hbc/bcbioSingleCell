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
#' @rdname diffExp
#' @name diffExp
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @inheritParams zingeR::zeroWeightsLS
#'
#' @param numerator Group of cells to use in the numerator of the contrast
#'   (e.g. treatment).
#' @param denominator Group of cells to use in the denominator of the contrast
#'   (e.g. control).
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
#' @return [DEGLRT] object.
#'
#' @examples
#' # seurat: expression in cluster 3 relative to cluster 2
#' ident2 <- WhichCells(pbmc_small, ident = 2L)
#' ident3 <- WhichCells(pbmc_small, ident = 3L)
#' lrt <- diffExp(
#'     pbmc_small,
#'     numerator = ident3,
#'     denominator = ident2,
#'     maxit = 100L)
NULL



# Constructors =================================================================
#' @importFrom basejump detectHPC initializeDirectory
#' @importFrom BiocParallel bpmapply
#' @importFrom edgeR calcNormFactors DGEList glmFit
#' @importFrom magrittr set_names
#' @importFrom Seurat WhichCells
#' @importFrom stats model.matrix
#' @importFrom zingeR glmWeightedF
.zingeR.edgeR.seurat <- function(  # nolint
    object,
    numerator,
    denominator,
    maxit = 1000L,
    dir = ".") {
    # Warn if user is running locally
    if (!detectHPC()) {
        inform(paste(
            "High-performance cluster environment not detected.",
            "It is strongly recommended to run differential expression",
            "as a job, since this can take hours to finish."
        ))
    }

    assert_is_all_of(object, "seurat")
    assert_is_character(numerator)
    assert_is_character(denominator)
    assert_are_disjoint_sets(numerator, denominator)
    dir <- initializeDirectory(dir)

    # Counts matrix
    cells <- c(numerator, denominator)
    counts <- counts(object, normalized = FALSE) %>%
        .[, cells]

    # Create a cell factor to define the group for `DGEList()`
    numeratorFactor <- replicate(
        n = length(numerator), expr = "numerator") %>%
        factor() %>%
        set_names(numerator)
    denominatorFactor <- replicate(
        n = length(denominator), expr = "denominator") %>%
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
    # levels(group)
    # table(group)
    # table(design[, 2L])

    # zingeR + edgeR analysis
    # Note that TMM needs to be consistently applied for both
    # `calcNormFactors()` and `zeroWeightsLS()`
    dge <- DGEList(counts, group = group)
    dge <- calcNormFactors(dge, method = "TMM")
    # dge[["samples"]]

    # This is the zingeR step that is computationally expensive
    weights <- zeroWeightsLS(
        counts = dge[["counts"]],
        design = design,
        maxit = maxit,
        normalization = "TMM")

    dge[["weights"]] <- weights
    dge <- estimateDisp(dge, design = design)

    # Plot biological coefficient of variation
    # plotBCV(dge)

    fit <- glmFit(dge, design = design)
    lrt <- glmWeightedF(fit, coef = 2L, independentFiltering = TRUE)

    # hist(lrt$table$PValue)
    # sum(lrt$table$padjFilter <= 0.05, na.rm = TRUE)

    lrt
}



# Methods ======================================================================
#' @rdname diffExp
#' @export
setMethod(
    "diffExp",
    signature("seurat"),
    .zingeR.edgeR.seurat)
