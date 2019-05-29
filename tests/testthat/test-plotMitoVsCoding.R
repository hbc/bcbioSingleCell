context("plotMitoVsCoding")

test_that("bcbioSingleCell", {
    x <- plotMitoVsCoding(indrops)
    expect_s3_class(x, "ggplot")
})
