context("subsetPerSample")

test_that("Unfiltered bcbioSingleCell", {
    expect_error(
        subsetPerSample(bcb),
        "`filterCells\\(\\)` hasn't been applied to this dataset"
    )
})

test_that("Filtered bcbioSingleCell", {
    expect_message(
        subsetPerSample(filtered, dir = "subsetPerSample"),
        "1 sample\\(s\\) matched: M1"
    )
    expect_identical(
        dir("subsetPerSample"),
        "M1.rda"
    )
    load("subsetPerSample/M1.rda")
    expect_identical(
        dim(M1),
        c(1000L, 243L)
    )
    unlink("subsetPerSample", recursive = TRUE)
})
