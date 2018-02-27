context("plotTopMarkers")

test_that("seurat", {
    top <- top[1L:2L, ,]
    expect_is(top, "grouped_df")
    x <- plotTopMarkers(seurat, topMarkers = top)
    expect_is(x, "list")
    expect_is(x[[1L]][[1L]], "ggplot")
})
