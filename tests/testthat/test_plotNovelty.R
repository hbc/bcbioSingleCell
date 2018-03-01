context("plotNovelty")

test_that("bcbioSingleCell", {
    p <- plotNovelty(bcb)
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotNovelty(seurat)
    expect_is(p, "ggplot")
})
