context("plotPCElbow")

test_that("plotPCElbow : seurat", {
    pcUse <- plotPCElbow(pbmc_small)
    expect_identical(pcUse, seq_len(14L))
})
