context("plotNovelty")

bcbFile <- system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell")
load(bcbFile)

test_that("plotNovelty", {
    p <- plotNovelty(bcb)
    expect_is(p, "ggplot")
})
