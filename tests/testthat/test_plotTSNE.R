context("plotTSNE")

test_that("seurat", {
    p <- plotTSNE(seurat)
    expect_is(p, "ggplot")
})
