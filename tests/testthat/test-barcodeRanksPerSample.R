context("barcodeRanksPerSample")

test_that("bcbioSingleCell", {
    skip_if_not(packageVersion("DropletUtils") >= "1.4")
    x <- barcodeRanksPerSample(indrops)
    expect_is(x, "list")
    expect_identical(
        names(x),
        "multiplexed_AAAAAAAA"
    )
    expect_s4_class(x[[1L]], "DataFrame")
    expect_identical(
        object = colnames(x[[1L]]),
        expected = c("rank", "total", "fitted")
    )
    ## DropletUtils currently returning numeric instead of integer.
    expect_equal(
        object = metadata(x[[1L]]),
        expected = list(
            knee = 2026L,
            inflection = 1205L
        )
    )
})
