context("plotPCA")

test_that("seurat", {
    p <- plotPCA(seurat)
    expect_is(p, "ggplot")
})
