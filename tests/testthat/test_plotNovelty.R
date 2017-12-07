context("plotNovelty")

test_that("plotNovelty", {
    p <- plotNovelty(bcb)
    expect_is(p, "ggplot")
})
