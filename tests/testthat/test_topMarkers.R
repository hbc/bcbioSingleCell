context("topMarkers")

test_that("grouped_df", {
    x <- topMarkers(all_markers_small)
    expect_is(x, "grouped_df")
    expect_identical(
        dplyr::group_vars(x),
        "cluster"
    )
    expect_identical(
        dim(x),
        c(50L, 18L)
    )
})
