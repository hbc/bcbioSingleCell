context("updateObject")

test_that("bcbioSingleCell", {
    x <- updateObject(indrops)
    expect_s4_class(x, "bcbioSingleCell")
})
