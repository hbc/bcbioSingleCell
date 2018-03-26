context("plotUMIsVsGenes")

test_that("bcbioSingleCell", {
    p <- plotUMIsVsGenes(bcb_small)
    expect_is(p, "ggplot")
})

test_that("seurat_small", {
    p <- plotUMIsVsGenes(seurat_small)
    expect_is(p, "ggplot")
})
