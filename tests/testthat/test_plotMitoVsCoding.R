context("plotMitoVsCoding")

load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))

test_that("bcbioSingleCell", {
    p <- plotMitoVsCoding(bcb)
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotMitoVsCoding(seurat)
    expect_is(p, "ggplot")
})

test_that("data.frame", {
    df <- metrics(bcb)
    expect_is(df, "data.frame")
    p <- plotMitoVsCoding(df)
    expect_is(p, "ggplot")
})
