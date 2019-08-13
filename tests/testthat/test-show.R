context("show")

test_that("bcbioSingleCell", {
    output <- capture.output(show(bcb))
    expect_true(grepl("^bcbioSingleCell", output[[1L]]))
})
