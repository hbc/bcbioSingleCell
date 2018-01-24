context("plotMarkerTSNE")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

symbol <- counts(seurat) %>%
    rownames() %>%
    .[[1L]]
ensgene <- bcbio(seurat, "gene2symbol") %>%
    .[which(.[["symbol"]] %in% symbol), "ensgene", drop = TRUE]

test_that("symbol", {
    expect_identical(symbol, "SCYL3")
    p <- plotMarkerTSNE(seurat, genes = symbol, format = "symbol")
    expect_is(p, "ggplot")
})

test_that("ensgene", {
    expect_identical(ensgene, "ENSG00000000457")
    p <- plotMarkerTSNE(seurat, genes = ensgene, format = "ensgene")
    expect_is(p, "ggplot")
})

test_that("data.frame", {
    df <- fetchTSNEExpressionData(seurat, genes = symbol)
    expect_is(df, "grouped_df")
    p <- plotMarkerTSNE(df)
    expect_is(p, "ggplot")
})
