context("plotReadsPerCell")

bcb <- examples[["bcb"]]

test_that("plotReadsPerCell", {
    p <- plotReadsPerCell(bcb)
    expect_is(p, "ggplot")
})
