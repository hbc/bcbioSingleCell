context("plotTopMarkers")

test_that("seurat", {
    topMarkers <- topMarkers[1L:2L, ,]
    expect_is(topMarkers, "grouped_df")
    plotlist <- plotTopMarkers(seurat, topMarkers)
    expect_is(plotlist, "list")
    expect_is(plotlist[[1L]][[1L]], "ggplot")
})
