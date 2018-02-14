context("detectOrganism")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("bcbioSingleCell", {
    expect_identical(
        detectOrganism(bcb),
        "Homo sapiens"
    )
})

test_that("seurat", {
    expect_identical(
        detectOrganism(seurat),
        "Homo sapiens"
    )
})
