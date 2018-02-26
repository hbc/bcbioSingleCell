context("plotMitoRatio")

test_that("bcbioSingleCell", {
    p <- plotMitoRatio(bcb)
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotMitoRatio(seurat)
    expect_is(p, "ggplot")
})
