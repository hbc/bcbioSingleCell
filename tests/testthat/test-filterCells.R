context("filterCells")

bcb <- calculateMetrics(bcb)

## FIXME now named "unknown"?
test_that("sampleNames", {
    expect_identical(
        object = sampleNames(bcb),
        expected = c("multiplexed_AAAAAAAA" = "rep_1")
    )
})

## Expecting an object with the same dimensions by default.
test_that("No filtering", {
    x <- filterCells(bcb)
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(dim(x), dim(bcb))
})

## Refer to the quality control R Markdown for actual recommended cutoffs.
## These are skewed, and designed to work with our minimal dataset.
test_that("Parameterized cutoff tests", {
    mapply(
        args = list(
            list("minCounts" = 2000L),
            list("maxCounts" = 2500L),
            list("minFeatures" = 45L),
            list("maxFeatures" = 49L),
            list("maxMitoRatio" = 0.1),
            list("minNovelty" = 0.5),
            list("minCellsPerFeature" = 95L)
        ),
        dim = list(
            c(50L, 35L),
            c(50L, 88L),
            c(50L, 95L),
            c(50L, 81L),
            c(50L, 22L),
            c(50L, 81L),
            c(45L, 100L)
        ),
        FUN = function(args, dim) {
            args[["object"]] <- bcb
            x <- do.call(what = filterCells, args = args)
            expect_s4_class(x, "bcbioSingleCell")
            expect_s4_class(metadata(x)[["filterCells"]], "SimpleList")
            expect_identical(metadata(x)[["subset"]], TRUE)
            expect_identical(dim(x), dim)
        },
        SIMPLIFY = FALSE
    )
})

test_that("Expected cutoff failure", {
    expect_error(
        object = filterCells(bcb, minCounts = Inf),
        expected = "No cells passed"
    )
})
