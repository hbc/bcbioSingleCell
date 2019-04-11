context("plotBarcodeRanks")

test_that("bcbioSingleCell", {
    x <- plotBarcodeRanks(indrops)
    expect_s3_class(x, "ggplot")
})
