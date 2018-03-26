context("sanitizeMarkers")

test_that("seurat_small", {
    x <- sanitizeMarkers(
        seurat_small,
        markers = seurat_all_markers_original)
    expect_is(x, "data.frame")
    annotable <- annotable(seurat_small)
    expect_identical(
        setdiff(colnames(x), colnames(annotable)),
        c("cluster",
          "pct1",
          "pct2",
          "avgLogFC",
          "pvalue",
          "padj")
    )
})
