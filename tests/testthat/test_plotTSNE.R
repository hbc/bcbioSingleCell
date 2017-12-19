context("plotTSNE")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("seurat", {
    p <- plotTSNE(seurat)
    expect_is(p, "ggplot")
})
