context("plotUMIsVsGenes")

test_that("bcbioSingleCell", {
    p <- plotUMIsVsGenes(bcb)
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotUMIsVsGenes(seurat)
    expect_is(p, "ggplot")
})
