context("plotNovelty")

bcb <- examples[["bcb"]]

test_that("plotNovelty", {
    p <- plotNovelty(bcb)
    expect_is(p, "ggplot")
})
