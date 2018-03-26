context("plotTSNE")

test_that("seurat_small", {
    p <- plotTSNE(seurat_small)
    expect_is(p, "ggplot")
})
