context("Data Functions")



# topBarcodes ==================================================================
test_that("topBarcodes", {
    # tibble
    x <- topBarcodes(indrops_small, return = "tibble")
    expect_identical(dplyr::group_vars(x), "sampleID")
    expect_identical(
        lapply(x, class),
        list(
            cellID = "character",
            sampleID = "factor",
            sampleName = "factor",
            nUMI = "integer"
        )
    )

    # list
    x <- topBarcodes(indrops_small, return = "list")
    expect_is(x, "list")
})
