context("plotQuantileHeatmap")

test_that("plotQuantileHeatmap : seurat", {
    p <- plotQuantileHeatmap(seurat_small)
    expect_is(p, "list")
    expect_identical(
        names(p),
        c("tree_row", "tree_col", "kmeans", "gtable")
    )
})
