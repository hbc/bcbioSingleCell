context("subsetPerSample")

test_that("subsetPerSample : bcbioSingleCell", {
    subsetPerSample(bcb_small, dir = "subsetPerSample")
    expect_identical(
        list.files("subsetPerSample"),
        "multiplexed_AAAAAAAA.rda"
    )
    load("subsetPerSample/multiplexed_AAAAAAAA.rda")
    expect_identical(
        dim(multiplexed_AAAAAAAA),
        c(500L, 500L)
    )
    unlink("subsetPerSample", recursive = TRUE)
})
