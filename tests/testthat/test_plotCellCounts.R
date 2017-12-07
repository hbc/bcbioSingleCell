context("plotCellCounts")

test_that("plotCellCounts", {
    p <- plotCellCounts(bcb)
    expect_is(p, "ggplot")
})
