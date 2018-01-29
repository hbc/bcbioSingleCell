context("plotMarkers")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "topMarkers.rda"),
    package = "bcbioSingleCell"))

test_that("symbol", {
    genes <- pull(topMarkers, "symbol") %>%
        .[[1L]]
    expect_identical(genes, "ACTC1")
    plotlist <- plotMarkers(seurat, genes = genes)
    expect_is(plotlist, "list")
    p <- plotlist[[1L]]
    expect_is(p, "ggplot")
})
