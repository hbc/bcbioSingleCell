context("plotReadsPerCell")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("bcbioSingleCell", {
    p <- plotReadsPerCell(bcb)
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotReadsPerCell(seurat)
    expect_is(p, "NULL")
})
