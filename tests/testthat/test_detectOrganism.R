context("detectOrganism")

test_that("bcbioSingleCell", {
    organism <- detectOrganism(bcb)
    expect_identical(
        organism,
        "Homo sapiens"
    )
})

test_that("seurat", {
    organism <- detectOrganism(seurat)
    expect_identical(
        organism,
        "Homo sapiens"
    )
})
