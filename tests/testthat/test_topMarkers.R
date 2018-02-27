context("topMarkers")

test_that("data.frame", {
    x <- topMarkers(seurat_all_markers)
    expect_is(x, "grouped_df")
    expect_identical(
        group_vars(x),
        "cluster"
    )
    expect_identical(
        dim(x),
        c(40L, 17L)
    )
})
