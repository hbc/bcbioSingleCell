context("plotMitoRatio")

bcb <- examples[["bcb"]]

test_that("plotMitoRatio", {
    p <- plotMitoRatio(bcb)
    expect_is(p, "ggplot")
})
