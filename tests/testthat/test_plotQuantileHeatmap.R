context("plotQuantileHeatmap")

test_that("plotQuantileHeatmap : seurat", {
    p <- plotQuantileHeatmap(seurat_small)
    # Check for pheatmap return
    expect_is(p, "list")
    expect_identical(
        names(p),
        c("quantiles", "plot")
    )
    expect_identical(
        names(p[["plot"]]),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
})
