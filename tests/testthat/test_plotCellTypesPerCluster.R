context("plotCellTypesPerCluster")

test_that("plotCellTypesPerCluster", {
    # Let's plot the first row, as an example
    x <- cellTypesPerCluster(known_markers_detected) %>%
        .[1L, , drop = FALSE]
    expect_identical(
        x[["symbol"]],
        "NCAM1"
    )
    p <- plotCellTypesPerCluster(
        seurat,
        cellTypesPerCluster = x)
    expect_is(p, "list")
    expect_is(p[[1L]][[1L]], "ggplot")
})
