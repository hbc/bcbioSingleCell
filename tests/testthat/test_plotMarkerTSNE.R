context("plotMarkerTSNE")

genes <- counts(seurat_small) %>%
    rownames() %>%
    .[[1L]]

test_that("seurat_small", {
    expect_identical(genes, "SCYL3")
    mean <- plotMarkerTSNE(seurat_small, genes = genes, expression = "mean")
    median <- plotMarkerTSNE(seurat_small, genes = genes, expression = "median")
    sum <- plotMarkerTSNE(seurat_small, genes = genes, expression = "sum")
    invisible(lapply(
        X = list(mean, median, sum),
        FUN = function(p) {
            expect_is(p, "ggplot")
        }
    ))
})

test_that("Invalid expression", {
    expect_error(
        plotMarkerTSNE(seurat_small, genes = genes, expression = "XXX"),
        "`expression` must contain: mean, median, sum"
    )
})

test_that("data.frame", {
    df <- fetchTSNEExpressionData(seurat_small, genes = genes)
    expect_is(df, "data.frame")
    p <- plotMarkerTSNE(df, genes = genes)
    expect_is(p, "ggplot")
})
