context("plotCellCounts")

test_that("bcbioSingleCell", {
    x <- plotCellCounts(indrops)
    expect_s3_class(x, "ggplot")
})
