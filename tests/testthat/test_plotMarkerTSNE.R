context("plotMarkerTSNE")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

genes <- counts(seurat) %>%
    rownames() %>%
    .[[1L]]

test_that("seurat", {
    expect_identical(genes, "SCYL3")
    p <- plotMarkerTSNE(seurat, genes = genes)
    expect_is(p, "ggplot")
})

test_that("data.frame", {
    df <- fetchTSNEExpressionData(seurat, genes = genes)
    expect_is(df, "data.frame")
    p <- plotMarkerTSNE(df, genes = genes)
    expect_is(p, "ggplot")
})
