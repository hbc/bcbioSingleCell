context("gene2symbol")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

colnames <- c("ensgene", "symbol")

test_that("bcbioSingleCell", {
    data <- gene2symbol(bcb)
    expect_is(data, "data.frame")
    expect_identical(colnames(data), colnames)
})

test_that("seurat", {
    data <- gene2symbol(seurat)
    expect_is(data, "data.frame")
    expect_identical(colnames(data), colnames)
})
