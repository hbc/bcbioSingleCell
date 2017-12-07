context("plotUMIsVsGenes")

test_that("plotUMIsVsGenes", {
    p <- plotUMIsVsGenes(bcb)
    expect_is(p, "ggplot")
})
