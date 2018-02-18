context("plotPCA")

load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))

test_that("seurat", {
    p <- plotPCA(seurat)
    expect_is(p, "ggplot")
})
