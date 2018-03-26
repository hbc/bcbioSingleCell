context("plotCellTypesPerCluster")

test_that("plotCellTypesPerCluster", {
    x <- cellTypesPerCluster(known_markers_small)[1L, , drop = FALSE]
    expect_identical(x[["geneName"]], "NCAM1")
    p <- plotCellTypesPerCluster(
        object = seurat_small,
        cellTypesPerCluster = x
    )
    expect_is(p, "list")
    expect_is(p[[1L]][[1L]], "ggplot")
})
