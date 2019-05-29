context("plotBarcodeRanks")

test_that("bcbioSingleCell", {
    skip_if_not(packageVersion("DropletUtils") >= "1.4")
    x <- plotBarcodeRanks(indrops)
    expect_s3_class(x, "ggplot")
})
