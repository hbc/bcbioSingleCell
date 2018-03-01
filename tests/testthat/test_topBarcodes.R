context("topBarcodes")

test_that("bcbioSingleCell", {
    data <- topBarcodes(bcb)
    expect_is(data, "data.frame")
    expect_identical(
        rownames(data)[[1L]],
        "run1_AGAGGATA_TGCTTCAT_GCAGGGTA"
    )
})

test_that("seurat", {
    data <- topBarcodes(seurat)
    expect_is(data, "data.frame")
    expect_identical(
        rownames(data)[[1L]],
        "M1_TGCTTCAT_GCAGGGTA"
    )
})
