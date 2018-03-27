context("plotMitoVsCoding")

test_that("bcbioSingleCell", {
    p <- plotMitoVsCoding(bcb_small)
    expect_is(p, "ggplot")
})

test_that("seurat_small", {
    p <- plotMitoVsCoding(seurat_small)
    expect_is(p, "ggplot")
})
