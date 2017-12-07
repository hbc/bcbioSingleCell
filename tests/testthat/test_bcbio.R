context("bcbio")

bcbFile <- system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell")
load(bcbFile)

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
