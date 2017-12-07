context("plotReadsPerCell")

test_that("plotReadsPerCell", {
    p <- plotReadsPerCell(bcb)
    expect_is(p, "ggplot")
})
