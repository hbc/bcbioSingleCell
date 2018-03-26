context("plotGenesPerCell")

test_that("bcbioSingleCell", {
    p <- plotGenesPerCell(bcb_small)
    expect_is(p, "ggplot")
})

test_that("seurat_small", {
    p <- plotGenesPerCell(seurat_small)
    expect_is(p, "ggplot")
})
