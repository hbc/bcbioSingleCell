context("plotTopMarkers")

load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
load(system.file("extdata/topMarkers.rda", package = "bcbioSingleCell"))

topMarkers <- topMarkers[1L:2L, ]

test_that("seurat", {
    expect_is(topMarkers, "grouped_df")
    plotlist <- plotTopMarkers(seurat, topMarkers)
    expect_is(plotlist, "list")
    expect_is(plotlist[[1L]][[1L]], "ggplot")
})
