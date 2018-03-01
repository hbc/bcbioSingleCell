context("plotUMIsPerCell")

test_that("bcbioSingleCell", {
    p <- plotUMIsPerCell(bcb)
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotUMIsPerCell(seurat)
    expect_is(p, "ggplot")
})
