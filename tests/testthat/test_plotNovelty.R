context("plotNovelty")

test_that("bcbioSingleCell", {
    p <- plotNovelty(bcb_small)
    expect_is(p, "ggplot")
})

test_that("seurat_small", {
    p <- plotNovelty(seurat_small)
    expect_is(p, "ggplot")
})
