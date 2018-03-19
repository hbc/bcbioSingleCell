context("aggregateReplicates")

test_that("aggregateReplicates", {
    x <- aggregateReplicates(bcb_small)
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(
        dim(x),
        c(500L, 286L)
    )
    map <- metadata(x)[["aggregateReplicates"]]
    expect_is(map, "factor")
    expect_identical(length(map), 500L)
    expect_identical(length(levels(map)), 286L)
})
