context("plotMarkers")

test_that("symbols", {
    genes <- pull(top, "geneName")[[1L]]
    expect_identical(genes, "ACTC1")
    plotlist <- plotMarkers(seurat_small, genes = genes)
    expect_is(plotlist, "list")
    p <- plotlist[[1L]]
    expect_is(p, "ggplot")
})
