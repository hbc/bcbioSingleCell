context("aggregateReplicates")

load(system.file(
    file.path("inst", "extdata", "bcb.rda"),
    package = "bcbioSingleCell"))

test_that("aggregateReplicates", {
    pooled <- suppressMessages(aggregateReplicates(bcb))
    expect_is(pooled, "bcbioSingleCell")
    expect_identical(
        dim(pooled),
        c(39297L, 10874L)
    )
    map <- metadata(pooled)[["aggregateReplicates"]]
    expect_is(map, "factor")
    expect_identical(
        length(map),
        18240L
    )
    expect_identical(
        length(levels(map)),
        10874L
    )
})
