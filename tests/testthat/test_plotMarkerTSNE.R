context("plotMarkerTSNE")

genes <- rownames(counts(pbmc_small))[[1L]]

test_that("Expression arguments", {
    args <- methodFormals("plotMarkerTSNE", "seurat") %>%
        .[["expression"]] %>%
        as.character() %>%
        .[-1L]
    invisible(lapply(args, function(arg) {
        p <- plotMarkerTSNE(pbmc_small, genes = genes, expression = arg)
        expect_is(p, "ggplot")
    }))
})
