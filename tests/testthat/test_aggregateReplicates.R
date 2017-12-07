context("aggregateReplicates")

bcbFile <- system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell")
load(bcbFile)

test_that("aggregateReplicates", {
    pooled <- suppressMessages(aggregateReplicates(bcb))
    expect_is(pooled, "bcbioSingleCell")
    expect_identical(
        dim(pooled),
        c(39297L, 36943L)
    )
    map <- metadata(pooled)[["aggregateReplicates"]]
    expect_is(map, "factor")
    expect_identical(
        length(map),
        62485L
    )
    expect_identical(
        length(levels(map)),
        36943L
    )
})
