context("plotMarkers")

test_that("symbol", {
    genes <- pull(top, "symbol")[[1L]]
    expect_identical(genes, "ACTC1")
    plotlist <- plotMarkers(seurat, genes = genes)
    expect_is(plotlist, "list")
    p <- plotlist[[1L]]
    expect_is(p, "ggplot")
})
