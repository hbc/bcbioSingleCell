context("topMarkers")

test_that("data.frame", {
    data <- topMarkers(seuratAllMarkers)
    expect_is(data, "grouped_df")
    expect_identical(
        group_vars(data),
        "cluster"
    )
    expect_identical(
        dim(data),
        c(40L, 17L)
    )
})
