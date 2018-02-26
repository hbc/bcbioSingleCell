# Differential Expression Using zingeR-edgeR
#
# Michael Steinbaugh
# 2018-02-26

library(BiocParallel)
library(edgeR)
library(zingeR)
library(Seurat)
library(bcbioSingleCell)

# ==============================================================================
# zingeR-edgeR
#
# We first perform a differential expression analysis using zingeR posterior
# probabilities and focussing on the count component of the ZINB model. We use
# edgeR to model the count component.
#
# The weights are estimated with the core function of zingeR, `zeroWeightsLS()`.
# It is important to be consistent with the normalization procedure, i.e. if you
# use TMM normalization for the analysis, you should also use it for estimating
# the zingeR posterior probabilities. The weights can then be provided in the
# DGEList object which will be used for dispersion estimation in the native
# edgeR `estimateDisp()` function.
#
# After estimation of the dispersions and zingeR posterior probabilities, the
# `glmWeightedF()` function is used for statistical inference. This is an
# adapted function from the `glmLRT()` function of edgeR. It uses an F-test for
# which the denominator degrees of freedom are by default adjusted according to
# the downweighting of excess zeros (`ZI = TRUE`). Also, independent filtering
# can be performed on the obtained p-values (`independentFiltering = TRUE`). We
# use the independent filtering strategy that was originally implemented in
# DESeq2. By default, the average fitted values are used as a filter criterion.
#
# https://goo.gl/4rTK1w
# ==============================================================================

data_dir <- path("dir", Sys.Date())
loadData(seurat, dir = data_dir)

# Get the cell cluster identities from the final Seurat object
ident <- slot(seurat, "ident") %>% levels()

dgelrt <- bpmapply(
    FUN = function(
        seurat,
        ident,
        envir = parent.frame(),
        dir = ".") {
        cells <- WhichCells(seurat, ident = ident)
        counts <- counts(seurat, normalized = FALSE)
        counts <- counts[, cells]

        # Get the primary interesting group, used to set up the contrast
        intgroup <- interestingGroups(seurat)[[1L]]
        group <- metrics(seurat) %>%
            .[cells, intgroup] %>%
            as.factor() %>%
            # Be sure to set the control as the reference level!
            relevel(ref = "control")

        # Check the counts
        table(group)

        # zingeR + edgeR analysis
        dge <- DGEList(counts, group = group)
        dge <- calcNormFactors(dge, method = "TMM")

        # Ensure the design matrix is set up correctly.
        # For comparing `mut` vs. `wt`, define `wt` as the reference level
        design <- model.matrix(~group)
        levels(group)  # wt, mut
        table(design[, 2L])

        # This is the zingeR step that is computationally expensive
        weights <- zeroWeightsLS(
            dge[["counts"]],
            design = design,
            # The `maxit` default of 200 is too low
            maxit = 1000L,
            normalization = "TMM")

        dge[["weights"]] <- weights
        dge <- estimateDisp(dge, design)

        plotBCV(dge)

        fit <- glmFit(dge, design)
        lrt <- glmWeightedF(fit, coef = 2L, independentFiltering = TRUE)

        # Save the individual DGELRT to disk, as a backup
        assignAndSaveData(
            name = paste("zinger_edger_lrt_cluster", ident, sep = "_"),
            object = lrt,
            envir = envir,
            dir = dir)

        lrt
    },
    ident = ident,
    MoreArgs = list(seurat = seurat, dir = data_dir),
    SIMPLIFY = FALSE,
    USE.NAMES = TRUE)
saveData(dgelrt, dir = dat_dir)
