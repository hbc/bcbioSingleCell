context("plotTopMarkers")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "topMarkers.rda"),
    package = "bcbioSingleCell"))
topMarkers <- topMarkers[1:2, ]

test_that("seurat", {
    expect_is(topMarkers, "grouped_df")
    plotlist <- plotTopMarkers(seurat, topMarkers)
    # FIXME This is returning empty lists
    expect_is(plotlist, "list")
    expect_is(plotlist[[1]][[1]], "ggplot")
})
