context("subsetPerSample")

test_that("bcbioSingleCell", {
    subsetPerSample(bcb_small, dir = "subsetPerSample")
    expect_identical(
        list.files("subsetPerSample"),
        c(
            "run1_AGAGGATA.rda",
            "run2_AGAGGATA.rda"
        )
    )
    load("subsetPerSample/run1_AGAGGATA.rda")
    expect_identical(
        dim(run1_AGAGGATA),
        c(500L, 286L)
    )
    unlink("subsetPerSample", recursive = TRUE)
})
