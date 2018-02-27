context("sanitizeMarkers")

test_that("seurat", {
    x <- sanitizeMarkers(
        seurat,
        markers = seurat_all_markers_original)
    expect_is(x, "data.frame")
    annotable <- annotable(seurat)
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
