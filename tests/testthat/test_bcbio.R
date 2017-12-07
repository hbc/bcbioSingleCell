context("bcbio")

bcb <- examples[["bcb"]]

test_that("bcbio", {
    slot <- bcbio(bcb)
    expect_is(slot, "SimpleList")
    expect_identical(
        lapply(slot, class),
        list(
            "cellularBarcodes" = "list"
        )
    )
})
