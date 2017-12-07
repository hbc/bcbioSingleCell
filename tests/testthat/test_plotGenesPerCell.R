context("plotGenesPerCell")

bcbFile <- system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell")
load(bcbFile)

test_that("plotGenesPerCell", {
    p <- plotGenesPerCell(bcb)
    expect_is(p, "ggplot")
})
