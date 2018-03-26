context("plotPCElbow")

test_that("seurat_small", {
    pcUse <- plotPCElbow(
        seurat_small,
        maxPct = 0.05,
        minCumPct = 0.8)
    expect_identical(
        pcUse,
        seq_len(10L)
    )
})
