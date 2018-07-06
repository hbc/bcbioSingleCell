# FIXME Removing zingeR support because it's not on Bioconductor yet



#' Differential Expression
#'
#' We now generally recommend the ZINB-WaVE method over zingeR, since it is
#' faster, and has been show to be more sensitive for most single-cell RNA-seq
#' datasets.
#'
#' @section zinbwave:
#' We are currently using an epsilon setting of `1e12`, as recommended by the
#' ZINB-WaVE integration paper.
#'
#' @section zingeR:
#' We first perform a differential expression analysis using zingeR posterior
#' probabilities and focussing on the count component of the ZINB model. The
#' weights are estimated with the core function of zingeR, `zeroWeightsLS()`. It
#' is important to be consistent with the normalization procedure, i.e. if you
#' use TMM normalization for the analysis, you should also use it for estimating
#' the zingeR posterior probabilities.
#'
#' @section edgeR:
#' After estimation of the dispersions and posterior probabilities, the
#' `glmWeightedF()` function is used for statistical inference. This is an
#' adapted function from the `glmLRT()` function of edgeR. It uses an F-test for
#' which the denominator degrees of freedom are by default adjusted according to
#' the downweighting of excess zeros (`ZI = TRUE`). Also, independent filtering
#' can be performed on the obtained p-values (`independentFiltering = TRUE`). We
#' use the independent filtering strategy that was originally implemented in
#' DESeq2. By default, the average fitted values are used as a filter criterion.
#'
#' @section DESeq2:
#' We're providing preliminary support for DESeq2 as the differential expression
#' caller. It is currently considerably slower for large datasets than edgeR.
#'
#' @seealso Consult the following sources for more information:
#' - [zinbwave paper](https://doi.org/10.1186/s13059-018-1406-4).
#' - [zinbwave vignette](https://bit.ly/2wtDdpS).
#' - [zinbwave-DESeq2 workflow](https://github.com/mikelove/zinbwave-deseq2).
#' - [zingeR vignette](https://goo.gl/4rTK1w).
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
#' @param caller Package to use for differential expression calling. Defaults to
#'   edgeR (faster for large datasets) but DESeq2 is also supported.
#' @param minCellsPerGene The minimum number of cells where a gene is expressed,
#' to pass low expression filtering. Set to `0` to disable (not recommended).
#' @param minCountsPerCell Minimum number of counts per cell for a gene to pass
#'   low expression filtering. The number of cells is defined by
#'   `minCellsPerGene`. Set to `0` to disable (not recommended).
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
#' @return
#' - `caller = "edgeR"`: `DEGLRT`.
#' - `caller = "DESeq2"`: Unshrunken `DESeqResults`. Use `lfcShrink()` if
#'   shrunken results are desired.
#'
#' @examples
#' # seurat ====
#' # Expression in cluster 3 relative to cluster 2
#' object <- seurat_small
#' numerator <- Seurat::WhichCells(object, ident = 3L)
#' glimpse(numerator)
#' denominator <- Seurat::WhichCells(object, ident = 2L)
#' glimpse(denominator)
#' x <- diffExp(
#'     object = object,
#'     numerator = numerator,
#'     denominator = denominator,
#'     minCellsPerGene = 5L,
#'     minCountsPerCell = 5L
#' )
NULL



.designFormula <- ~group



.zinbwave <- function(object) {
    message("Running zinbwave")
    stopifnot(is(object, "SingleCellExperiment"))
    object <- as(object, "SingleCellExperiment")
    # zinbFit doesn't support `dgCMatrix``, so coerce counts to matrix
    assays(object) <- list("counts" = as.matrix(counts(object)))
    print(system.time({
        zinb <- zinbwave::zinbwave(
            Y = object,
            K = 0L,
            BPPARAM = SerialParam(),
            epsilon = 1e12
        )
    }))
    stopifnot(is(zinb, "SingleCellExperiment"))
    assert_is_factor(zinb[["group"]])
    assert_is_matrix(metadata(zinb)[["design"]])
    zinb
}



# Van De Berge and Perraudeau and others have shown the LRT may perform better
# for null hypothesis testing, so we use the LRT. In order to use the Wald test,
# it is recommended to set `useT = TRUE`.
#
# For UMI data, for which the expected counts may be very low, the likelihood
# ratio test implemented in nbinomLRT should be used.
#
# DESeq2 supports `weights` in assays automatically.
.zinbwave.DESeq2 <- function(object) {  # nolint
    stopifnot(is(object, "SingleCellExperiment"))
    zinb <- .zinbwave(object)
    # DESeq2 ===================================================================
    message("Running DESeq2")
    print(system.time({
        dds <- DESeq2::DESeqDataSet(se = zinb, design = .designFormula)
        dds <- DESeq2::DESeq(
            object = dds,
            test = "LRT",
            reduced = ~ 1L,
            sfType = "poscounts",
            minmu = 1e-6,
            minReplicatesForReplace = Inf
        )
        # We already performed low count filtering
        res <- DESeq2::results(dds, independentFiltering = FALSE)
    }))
    res
}



.zinbwave.edgeR <- function(object) {  # nolint
    stopifnot(is(object, "SingleCellExperiment"))
    zinb <- .zinbwave(object)
    # edgeR ====================================================================
    message("Running edgeR")
    counts <- as.matrix(counts(zinb))
    weights <- assay(zinb, "weights")
    assert_is_matrix(weights)
    design <- metadata(object)[["design"]]
    assert_is_matrix(design)
    group <- object[["group"]]
    assert_is_factor(group)
    print(system.time({
        dge <- edgeR::DGEList(counts, group = group)
        dge <- edgeR::calcNormFactors(dge)
        dge[["weights"]] <- weights
        dge <- edgeR::estimateDisp(dge, design = design)
        fit <- edgeR::glmFit(dge, design = design)
        # We already performed low count filtering
        lrt <- zinbwave::glmWeightedF(
            glmfit = fit,
            coef = 2L,
            independentFiltering = FALSE
        )
    }))
    lrt
}



# Note that TMM needs to be consistently applied for both
# `calcNormFactors()` and `zeroWeightsLS()`.
.zingeR.edgeR <- function(  # nolint
    object,
    maxit = 1000L
) {
    stopifnot(is(object, "SingleCellExperiment"))
    counts <- as.matrix(counts(object))
    design <- metadata(object)[["design"]]
    assert_is_matrix(design)
    # zingeR ===================================================================
    message("Running zingeR")
    print(system.time({
        weights <- zingeR::zeroWeightsLS(
            counts = counts,
            design = design,
            maxit = maxit,
            normalization = "TMM"
        )
    }))
    # edgeR ====================================================================
    message("Running edgeR")
    group <- object[["group"]]
    assert_is_factor(group)
    print(system.time({
        dge <- edgeR::DGEList(counts, group = group)
        dge[["weights"]] <- weights
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        dge <- edgeR::estimateDisp(dge, design = design)
        fit <- edgeR::glmFit(dge, design = design)
        # We already performed low count filtering
        lrt <- zingeR::glmWeightedF(
            glmfit = fit,
            coef = 2L,
            independentFiltering = FALSE
        )
    }))
    lrt
}



.zingeR.DESeq2 <- function(  # nolint
    object,
    maxit = 1000L
) {
    stopifnot(is(object, "SingleCellExperiment"))
    counts <- as.matrix(counts(object))
    se <- as(object, "RangedSummarizedExperiment")
    assays(se) <- list(counts = counts)
    # zingeR ===================================================================
    message("Running zingeR")
    design <- metadata(object)[["design"]]
    assert_is_matrix(design)
    print(system.time({
        weights <- zingeR::zeroWeightsLS(
            counts = counts,
            design = design,
            maxit = maxit,
            normalization = "DESeq2_poscounts",
            colData = colData(object),
            designFormula = .designFormula
        )
    }))
    # DESeq2 ===================================================================
    message("Running DESeq2")
    print(system.time({
        dds <- DESeq2::DESeqDataSet(se = se, design = .designFormula)
        assays(dds)[["weights"]] <- weights
        dds <- DESeq2::estimateSizeFactors(dds, type = "poscounts")
        dds <- DESeq2::estimateDispersions(dds)
        # LRT is recommended over Wald for UMI counts.
        # May want to switch to that approach here.
        # See code in zinbwave method.
        dds <- DESeq2::nbinomWaldTest(
            object = dds,
            betaPrior = TRUE,
            useT = TRUE,
            df = rowSums(weights) - 2L
        )
        # We already performed low count filtering
        res <- DESeq2::results(dds, independentFiltering = FALSE)
    }))
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
        zeroWeights = c("zinbwave", "zingeR"),
        caller = c("edgeR", "DESeq2"),
        minCellsPerGene = 25L,
        minCountsPerCell = 5L
    ) {
        object <- as(object, "SingleCellExperiment")
        assert_is_character(numerator)
        assert_is_character(denominator)
        assert_are_disjoint_sets(numerator, denominator)
        zeroWeights <- match.arg(zeroWeights)
        if (zeroWeights == "zinbwave") {
            requireNamespace(
                package = "zinbwave",
                versionCheck = list(
                    op = ">=",
                    version = package_version("1.0")
                )
            )
        } else if (zeroWeights == "zingeR") {
            requireNamespace(
                package = "zingeR",
                versionCheck = list(
                    op = ">=",
                    version = package_version("0.1")
                )
            )
        }
        caller <- match.arg(caller)
        if (caller == "DESeq2") {
            requireNamespace(
                package = "DESeq2",
                versionCheck = list(
                    op = ">=",
                    version = package_version("1.20")
                )
            )
        } else if (zeroWeights == "edgeR") {
            requireNamespace(
                package = "edgeR",
                versionCheck = list(
                    op = ">=",
                    version = package_version("3.22")
                )
            )
        }
        assertIsAnImplicitInteger(minCountsPerCell)
        assertIsAnImplicitInteger(minCellsPerGene)

        message(paste(
            "Performing differential expression with",
            paste(zeroWeights, caller, sep = "-")
        ))

        # Subset the SCE object to contain the desired cells
        cells <- c(numerator, denominator)
        object <- object[, cells]

        message(paste(
            "Applying low gene count filter",
            paste(
                "Requiring at least",
                minCellsPerGene,
                "cells with counts of",
                minCountsPerCell,
                "or more per gene"
            ),
            sep = "\n"
        ))
        genes <- rowSums(counts(object) >= minCountsPerCell) >= minCellsPerGene
        # Early return NULL if no genes pass
        if (!any(genes)) {
            warning("No genes passed the low count filter")
            return(NULL)
        }
        object <- object[genes, ]

        # Create a cell factor to define the group for `DGEList()`
        numeratorFactor <- replicate(
            n = length(numerator),
            expr = "numerator"
        ) %>%
            as.factor() %>%
            set_names(numerator)
        denominatorFactor <- replicate(
            n = length(denominator),
            expr = "denominator"
        ) %>%
            as.factor() %>%
            set_names(denominator)
        group <- factor(c(
            as.character(numeratorFactor),
            as.character(denominatorFactor)
        ))
        names(group) <- c(names(numeratorFactor), names(denominatorFactor))
        # Ensure denominator is set as reference
        group <- relevel(group, ref = "denominator")
        object[["group"]] <- group

        # Set up the design matrix
        design <- model.matrix(~group)
        metadata(object)[["design"]] <- design

        fun <- get(paste("", zeroWeights, caller, sep = "."))
        fun(object)
    }
)



#' @rdname diffExp
#' @export
setMethod(
    "diffExp",
    signature("seurat"),
    getMethod("diffExp", "SingleCellExperiment")
)
