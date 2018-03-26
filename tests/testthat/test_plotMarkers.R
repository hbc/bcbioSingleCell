context("plotMarkers")

test_that("symbols", {
    genes <- head(rownames(seurat_small), n = 2L)
    plotlist <- plotMarkers(
        object = seurat_small,
        genes = genes
    )
    expect_is(plotlist, "list")
    expect_identical(length(plotlist), 2L)
    expect_identical(names(plotlist), names(genes))
    p <- plotlist[[1L]]
    expect_is(p, "ggplot")
})
