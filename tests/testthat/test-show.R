context("show")

test_that("bcbioSingleCell", {
    output <- capture.output(show(indrops))
    expect_true(grepl("^bcbioSingleCell", output[[1L]]))
})
