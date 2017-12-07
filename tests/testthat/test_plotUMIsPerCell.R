context("plotUMIsPerCell")

bcb <- examples[["bcb"]]

test_that("plotUMIsPerCell", {
    p <- plotUMIsPerCell(bcb)
    expect_is(p, "ggplot")
})
