context("plotMitoRatio")

bcbFile <- system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell")
load(bcbFile)

test_that("plotMitoRatio", {
    p <- plotMitoRatio(bcb)
    expect_is(p, "ggplot")
})
