context("plotCellCounts")

test_that("bcbioSingleCell", {
    p <- plotCellCounts(bcb_small)
    expect_is(p, "ggplot")
})

test_that("seurat_small", {
    p <- plotCellCounts(seurat_small)
    expect_is(p, "ggplot")
})
