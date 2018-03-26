context("plotZerosVsDepth")

test_that("bcbioSingleCell", {
    p <- plotZerosVsDepth(bcb_small)
    expect_is(p, "ggplot")
})

test_that("dgCMatrix", {
    counts <- counts(bcb_small)
    metrics <- metrics(bcb_small)
    p <- plotZerosVsDepth(counts, metrics = metrics)
    expect_is(p, "ggplot")
})
