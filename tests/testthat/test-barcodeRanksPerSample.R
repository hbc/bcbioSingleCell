context("barcodeRanksPerSample")

test_that("bcbioSingleCell", {
    x <- barcodeRanksPerSample(indrops)
    expect_is(x, "list")
    expect_identical(
        names(x),
        "multiplexed_AAAAAAAA"
    )
    expect_identical(
        object = names(x[[1L]]),
        expected = c("rank", "total", "fitted", "knee", "inflection")
    )
    # DropletUtils currently returning numeric instead of integer for these.
    expect_equal(
        object = x[["multiplexed_AAAAAAAA"]][c("knee", "inflection")],
        expected = list(
            knee = c(CTCTATAG_GCAAAGCC = 2026L),
            inflection = c(AATCGAAG_CCCAAGCA = 1205L)
        )
    )
})
