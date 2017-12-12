context("plotMarkers")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "topMarkers.rda"),
    package = "bcbioSingleCell"))

test_that("ensgene", {
    ensgene <- topMarkers$ensgene[[1]]
    expect_identical(ensgene, "ENSG00000159251")
    plotlist <- plotMarkers(seurat, genes = ensgene, format = "ensgene")
    expect_is(plotlist, "list")
    expect_is(plotlist[[1L]], "ggplot")
})

test_that("symbol", {
    symbol <- topMarkers$symbol[[1]]
    expect_identical(symbol, "ACTC1")
    p <- plotMarkers(seurat, genes = symbol, format = "symbol")
    expect_is(plotlist, "list")
    expect_is(plotlist[[1L]], "ggplot")
})
