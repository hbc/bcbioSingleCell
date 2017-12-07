context("plotCellCounts")

bcbFile <- system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell")
load(bcbFile)

test_that("plotCellCounts", {
    p <- plotCellCounts(bcb)
    expect_is(p, "ggplot")
})
