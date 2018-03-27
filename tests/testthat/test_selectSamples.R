context("selectSamples")

test_that("selectSamples : bcbioSingleCell", {
    x <- selectSamples(bcb_small, sampleName = "rep_1")
    expect_s4_class(x, "bcbioSingleCell")
    expect_true(metadata(x)[["selectSamples"]])
    expect_identical(
        dim(x),
        c(500L, 500L)
    )
    expect_identical(
        sampleData(x)[["sampleID"]],
        factor("multiplexed_AAAAAAAA")
    )
    # Ensure that bcbio cellular barcodes get updated correctly
    cb <- metadata(x)[["cellularBarcodes"]]
    expect_is(cb, "list")
    ids1 <- sampleData(x) %>%
        .[["sampleID"]] %>%
        as.character()
    ids2 <- names(cb)
    expect_identical(ids1, ids2)
    # Check that tibbles are stored inside list
    expect_identical(
        lapply(cb, class) %>%
            unlist() %>%
            unique(),
        c("tbl_df", "tbl", "data.frame")
    )
})

test_that("selectSamples : Match failure", {
    expect_error(
        selectSamples(bcb_small, sampleID = "XXX"),
        "sampleID metadata column doesn't contain XXX"
    )
})
