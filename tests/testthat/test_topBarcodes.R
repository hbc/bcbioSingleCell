context("topBarcodes")

test_that("bcbioSingleCell", {
    x <- topBarcodes(bcb_small)
    expect_is(x, "data.frame")
    expect_identical(
        rownames(x)[[1L]],
        "run1_AGAGGATA_TGCTTCAT_GCAGGGTA"
    )
})

test_that("seurat_small", {
    x <- topBarcodes(pbmc_small)
    expect_is(x, "data.frame")
    expect_identical(
        rownames(x)[[1L]],
        "GACATTCTCCACCT"
    )
})
