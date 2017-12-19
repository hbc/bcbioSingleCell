context("plotQC")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("bcbioSingleCell", {
    # Example dataset doesn't have a cellular barcode cutoff because we removed
    # the bcbio commands log file (which conflicts with Travis CI)
    p <- suppressWarnings(plotQC(bcb))
    expect_is(p, "ggplot")
})

test_that("seurat", {
    p <- plotQC(seurat)
    expect_is(p, "ggplot")
})
