context("plotUMIsVsGenes")

bcbFile <- system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell")
load(bcbFile)

test_that("plotUMIsVsGenes", {
    p <- plotUMIsVsGenes(bcb)
    expect_is(p, "ggplot")
})
