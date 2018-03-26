context("plotUMIsPerCell")

test_that("bcbioSingleCell", {
    p <- plotUMIsPerCell(bcb_small)
    expect_is(p, "ggplot")
})

test_that("seurat_small", {
    p <- plotUMIsPerCell(seurat_small)
    expect_is(p, "ggplot")
})
