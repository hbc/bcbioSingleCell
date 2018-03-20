context("topMarkers")

test_that("grouped_df", {
    x <- topMarkers(all_markers_small)
    expect_is(x, "grouped_df")
    expect_identical(
        group_vars(x),
        "cluster"
    )
    expect_identical(
        dim(x),
        c(40L, 18L)
    )
})
