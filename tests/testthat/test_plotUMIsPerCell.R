context("plotUMIsPerCell")

test_that("plotUMIsPerCell : bcbioSingleCell", {
    p <- plotUMIsPerCell(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotUMIsPerCell : seurat", {
    p <- plotUMIsPerCell(seurat_small)
    expect_is(p, "ggplot")
})
