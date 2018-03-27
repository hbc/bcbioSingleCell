context("plotTSNE")

test_that("plotTSNE : seurat", {
    p <- plotTSNE(pbmc_small)
    expect_is(p, "ggplot")
})
