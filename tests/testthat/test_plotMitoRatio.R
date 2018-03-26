context("plotMitoRatio")

test_that("bcbioSingleCell", {
    p <- plotMitoRatio(bcb_small)
    expect_is(p, "ggplot")
})

test_that("seurat_small", {
    p <- plotMitoRatio(seurat_small)
    expect_is(p, "ggplot")
})
