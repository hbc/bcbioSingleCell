context("plotPCElbow")

load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))

test_that("seurat", {
    pcUse <- plotPCElbow(
        seurat,
        maxPct = 0.05,
        minCumPct = 0.8)
    expect_identical(
        pcUse,
        seq_len(10L)
    )
})
