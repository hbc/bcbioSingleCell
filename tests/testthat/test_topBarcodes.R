context("topBarcodes")

test_that("topBarcodes : bcbioSingleCell", {
    x <- topBarcodes(bcb_small)
    expect_is(x, "data.frame")
    expect_identical(
        rownames(x)[[1L]],
        "multiplexed_AAAAAAAA_CTAGCACG_AATCGGGT"
    )
})

test_that("topBarcodes : seurat_small", {
    x <- topBarcodes(pbmc_small)
    expect_is(x, "data.frame")
    expect_identical(
        rownames(x)[[1L]],
        "GACATTCTCCACCT"
    )
})
