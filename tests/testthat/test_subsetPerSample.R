context("subsetPerSample")

load(system.file(
    file.path("extdata", "filtered.rda"),
    package = "bcbioSingleCell"))

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
