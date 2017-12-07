context("plotReadsPerCell")

bcbFile <- system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell")
load(bcbFile)

test_that("plotReadsPerCell", {
    p <- plotReadsPerCell(bcb)
    expect_is(p, "ggplot")
})
