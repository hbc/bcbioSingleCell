context("plotUMIsPerCell")

bcbFile <- system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell")
load(bcbFile)

test_that("plotUMIsPerCell", {
    p <- plotUMIsPerCell(bcb)
    expect_is(p, "ggplot")
})
