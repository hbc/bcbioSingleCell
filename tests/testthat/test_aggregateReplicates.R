context("aggregateReplicates")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))

test_that("aggregateReplicates", {
    pooled <- suppressMessages(aggregateReplicates(bcb))
    expect_is(pooled, "bcbioSingleCell")
    expect_identical(
        dim(pooled),
        c(1000L, 286L)
    )
    map <- metadata(pooled)[["aggregateReplicates"]]
    expect_is(map, "factor")
    expect_identical(length(map), 500L)
    expect_identical(length(levels(map)), 286L)
})
