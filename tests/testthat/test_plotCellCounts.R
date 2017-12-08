context("plotCellCounts")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("bcbioSingleCell", {
    p <- plotCellCounts(bcb)
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotCellCounts(seurat)
    expect_is(p, "ggplot")
})
