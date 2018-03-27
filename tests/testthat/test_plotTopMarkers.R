context("plotTopMarkers")

test_that("plotTopMarkers : seurat", {
    top <- topMarkers(all_markers_small, n = 1L)
    expect_is(top, "grouped_df")
    expect_identical(nrow(top), 5L)
    top <- top[seq_len(2L), ]
    x <- plotTopMarkers(seurat_small, topMarkers = top)
    expect_is(x, "list")
    expect_is(x[[1L]][[1L]], "ggplot")
})
