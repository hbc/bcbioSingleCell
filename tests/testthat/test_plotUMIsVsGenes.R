context("plotUMIsVsGenes")

test_that("plotUMIsVsGenes : bcbioSingleCell", {
    p <- plotUMIsVsGenes(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotUMIsVsGenes : seurat", {
    p <- plotUMIsVsGenes(seurat_small)
    expect_is(p, "ggplot")

    p <- plotUMIsVsGenes(pbmc_small)
    expect_is(p, "ggplot")
})
