context("plotPCA")

test_that("seurat_small", {
    p <- plotPCA(seurat_small)
    expect_is(p, "ggplot")
})
