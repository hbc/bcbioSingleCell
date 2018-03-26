context("plotMarkerTSNE")

genes <- rownames(counts(pbmc_small))[[1L]]

test_that("Expression arguments", {
    args <- methodFormals("plotMarkerTSNE", "seurat") %>%
        .[["expression"]] %>%
        as.character() %>%
        .[-1L]
    invisible(lapply(args, function(arg) {
        p <- plotMarkerTSNE(pbmc_small, genes = genes, expression = arg)
        expect_is(p, "ggplot")
    }))
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
