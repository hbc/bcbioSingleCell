context("plotQC")

bcbFile <- system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell")
load(bcbFile)

test_that("plotQC", {
    p <- plotQC(bcb, return = "grid")
    expect_is(p, "ggplot")
    list <- plotQC(bcb, return = "list")
    expect_is(list, "list")
    md <- capture.output(plotQC(bcb, return = "markdown"))
    expect_is(md, "character")
})
