context("plotPCA")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("seurat", {
    p <- plotPCA(seurat)
    expect_is(p, "ggplot")
})
