context("plotGenesPerCell")

bcb <- examples[["bcb"]]

test_that("plotGenesPerCell", {
    p <- plotGenesPerCell(bcb)
    expect_is(p, "ggplot")
})
