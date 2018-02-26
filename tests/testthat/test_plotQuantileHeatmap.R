context("plotQuantileHeatmap")

test_that("seurat", {
    p <- plotQuantileHeatmap(seurat)
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
