context("plotPCElbow")

test_that("plotPCElbow : seurat", {
    x <- plotPCElbow(seurat_small)
    expect_identical(x, seq_len(9L))

    x <- plotPCElbow(pbmc_small)
    expect_identical(x, seq_len(11L))
})
