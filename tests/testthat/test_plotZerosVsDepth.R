context("plotZerosVsDepth")

test_that("bcbioSingleCell", {
    p <- plotZerosVsDepth(bcb)
    expect_is(p, "ggplot")
})

test_that("dgCMatrix", {
    counts <- counts(bcb)
    metrics <- metrics(bcb)
    p <- plotZerosVsDepth(counts, metrics = metrics)
    expect_is(p, "ggplot")
})
