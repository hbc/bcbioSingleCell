context("plotGenesPerCell")

test_that("bcbioSingleCell", {
    p <- plotGenesPerCell(bcb)
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotGenesPerCell(seurat)
    expect_is(p, "ggplot")
})
