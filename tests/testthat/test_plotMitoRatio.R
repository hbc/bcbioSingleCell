context("plotMitoRatio")

load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))

test_that("bcbioSingleCell", {
    p <- plotMitoRatio(bcb)
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotMitoRatio(seurat)
    expect_is(p, "ggplot")
})
