context("detectOrganism")

load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))

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
