context("plotGenesPerCell")

test_that("plotGenesPerCell", {
    p <- plotGenesPerCell(bcb)
    expect_is(p, "ggplot")
})
