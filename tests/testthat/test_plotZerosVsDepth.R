context("plotZerosVsDepth")

test_that("plotZerosVsDepth : bcbioSingleCell", {
    p <- plotZerosVsDepth(bcb_small)
    expect_is(p, "ggplot")
})

test_that("plotZerosVsDepth : seurat", {
    p <- plotZerosVsDepth(pbmc_small)
    expect_is(p, "ggplot")
})
