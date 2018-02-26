context("plotCellCounts")

test_that("bcbioSingleCell", {
    p <- plotCellCounts(bcb)
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotCellCounts(seurat)
    expect_is(p, "ggplot")
})
