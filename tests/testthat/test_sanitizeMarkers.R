context("sanitizeMarkers")

test_that("seurat", {
    data <- sanitizeMarkers(
        seurat,
        markers = seuratAllMarkersOriginal)
    expect_is(data, "data.frame")
    annotable <- annotable(seurat)
    expect_identical(
        setdiff(colnames(data), colnames(annotable)),
        c("cluster",
          "pct1",
          "pct2",
          "avgLogFC",
          "pvalue",
          "padj")
    )
})
