context("plotCellTypesPerCluster")

load(system.file(
    file.path("extdata", "knownMarkersDetected.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("plotCellTypesPerCluster", {
    # Let's plot the first row, as an example
    cellTypesPerCluster <- cellTypesPerCluster(knownMarkersDetected) %>%
        .[1L, , drop = FALSE]
    expect_identical(
        cellTypesPerCluster[["symbol"]],
        "NCAM1"
    )
    plotlist <- plotCellTypesPerCluster(
        seurat,
        cellTypesPerCluster = cellTypesPerCluster)
    expect_is(plotlist, "list")
    expect_is(plotlist[[1L]][[1L]], "ggplot")
})
