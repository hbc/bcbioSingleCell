context("plotMarkerTSNE")

load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))

genes <- counts(seurat) %>%
    rownames() %>%
    .[[1L]]

test_that("seurat", {
    expect_identical(genes, "SCYL3")
    mean <- plotMarkerTSNE(seurat, genes = genes, expression = "mean")
    median <- plotMarkerTSNE(seurat, genes = genes, expression = "median")
    sum <- plotMarkerTSNE(seurat, genes = genes, expression = "sum")
    invisible(lapply(
        X = list(mean, median, sum),
        FUN = function(p) {
            expect_is(p, "ggplot")
        }
    ))
})

test_that("Invalid expression", {
    expect_error(
        plotMarkerTSNE(seurat, genes = genes, expression = "XXX"),
        "`expression` must contain: mean, median, sum"
    )
})

test_that("data.frame", {
    df <- fetchTSNEExpressionData(seurat, genes = genes)
    expect_is(df, "data.frame")
    p <- plotMarkerTSNE(df, genes = genes)
    expect_is(p, "ggplot")
})
