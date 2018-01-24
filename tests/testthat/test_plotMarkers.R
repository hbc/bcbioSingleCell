context("plotMarkers")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "topMarkers.rda"),
    package = "bcbioSingleCell"))

test_that("ensgene", {
    ensgene <- topMarkers[["ensgene"]][[1L]]
    expect_identical(ensgene, "ENSG00000159251")
    plotlist <- plotMarkers(seurat, genes = ensgene, format = "ensgene")
    expect_is(plotlist, "list")
    p <- plotlist[[1L]]
    expect_is(p, "ggplot")
})

test_that("symbol", {
    symbol <- topMarkers[["symbol"]][[1L]]
    expect_identical(symbol, "ACTC1")
    plotlist <- plotMarkers(seurat, genes = symbol, format = "symbol")
    expect_is(plotlist, "list")
    p <- plotlist[[1L]]
    expect_is(p, "ggplot")
})
