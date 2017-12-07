context("plotMitoRatio")

test_that("plotMitoRatio", {
    p <- plotMitoRatio(bcb)
    expect_is(p, "ggplot")
})
