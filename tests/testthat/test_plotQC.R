context("plotQC")

load(system.file(
    file.path("inst", "extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("inst", "extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("bcbioSingleCell", {
    p <- plotQC(bcb, return = "grid")
    expect_is(p, "ggplot")
    list <- plotQC(bcb, return = "list")
    expect_is(list, "list")
    md <- capture.output(plotQC(bcb, return = "markdown"))
    expect_is(md, "character")
})

test_that("seurat", {
    p <- plotQC(seurat, return = "grid")
    expect_is(p, "ggplot")
})
