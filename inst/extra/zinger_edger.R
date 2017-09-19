# Differential expression using zingeR + edgeR
# Michael Steinbaugh
# 2017-09-19

library(pbapply)
library(bcbioSingleCell)
library(Seurat)
library(zingeR)

# ==============================================================================
# zingeR-edgeR
#
# We first perform a differential expression analysis using zingeR posterior
# probabilities and focussing on the count component of the ZINB model. We use
# edgeR to model the count component.
#
# The weights are estimated with the core function of zingeR, zeroWeightsLS. It
# is important to be consistent with the normalization procedure, i.e. if you
# use TMM normalization for the analysis, you should also use it for estimating
# the zingeR posterior probabilities. In a next section, we demonstrate a
# zero-inflated NB analysis using DESeq2 where we re-estimate the weights using
# the proper DESeq2 normalization. The user can also input normalization factors
# from any normalization procedure of interest by using the normFactors argument
# in zeroWeightsLS.
#
# The weights can then be provided in the DGEList object which will be used for
# dispersion estimation in the native edgeR estimateDisp function.
#
# After estimation of the dispersions and zingeR posterior probabilities, the
# glmWeightedF function is used for statistical inference. This is an adapted
# function from the glmLRT function of edgeR. It uses an F-test for which the
# denominator degrees of freedom are by default adjusted according to the
# downweighting of excess zeros (argument ZI=TRUE). Also, independent filtering
# can be performed on the obtained p-values (argument
# independentFiltering=TRUE). We use the independent filtering strategy that was
# originally implemented in DESeq2. By default, the average fitted values are
# used as a filter criterion.
# 
# https://github.com/statOmics/zingeR/blob/master/vignettes/zingeRVignette_v2.Rmd
# ==============================================================================

loadData(bcb, seurat)
intgroup <- interestingGroups(bcb)[[1]]

# Get the cell cluster identities from the final Seurat object
ident <- levels(seurat@ident)

zinger_edger_lrt <- pblapply(seq_along(ident), function(a) {
    cells <- WhichCells(seurat, ident = ident[[a]])
    length(cells)
    
    bcb_subset <- bcb[, cells]
    counts <- assay(bcb_subset)
    
    group <- metrics(bcb_subset) %>%
        .[, intgroup] %>%
        as.factor %>%
        # Be sure to set wild-type (`wt`) as the reference level
        relevel(ref = "wt")
    
    # Check the counts
    table(group)
    
    # zingeR + edgeR analysis
    dge <- DGEList(counts, group = group)
    dge <- calcNormFactors(dge, method = "TMM")
    
    # Ensure the design matrix is set up correctly.
    # For comparing `mut` vs. `wt`, define `wt` as the reference level
    design <- model.matrix(~group)
    levels(group)  # c("wt", "mut")
    table(design[, 2])
    
    # This is the zingeR step that is computationally expensive
    weights <- zeroWeightsLS(
        dge$counts,
        design = design,
        normalization = "TMM")
    
    dge$weights <- weights
    dge <- estimateDisp(dge, design)

    plotBCV(dge)
    
    fit <- glmFit(dge, design)
    lrt <- glmWeightedF(fit, coef = 2, independentFiltering = TRUE)

    # hist(lrt$table$PValue)
    # sum(lrt$table$padjFilter <= 0.05, na.rm = TRUE)

    # DEG data.frame can be returned with:
    # lrt$table

    # Save the DGELRT to disk
    name <- paste0("zinger_edger_lrt_cluster_", ident[[a]])
    assignAndSaveData(name, lrt)
    
    lrt
})

# Optionally save all DGELRT in a single list
saveData(zinger_edger_lrt)
