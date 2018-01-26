context("interestingGroups")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("bcbioSingleCell", {
    x <- interestingGroups(bcb)
    expect_identical(x, "sampleName")
})

test_that("Assignment of undefined interesting group", {
    error <- "Interesting groups not defined in metadata: XXX"
    expect_error(
        interestingGroups(bcb) <- "XXX",
        error
    )
    expect_error(
        interestingGroups(seurat) <- "XXX",
        error
    )
})

test_that("seurat", {
    x <- interestingGroups(seurat)
    expect_identical(x, "sampleName")
})

test_that("seurat without bcbio metadata", {
    slot(seurat, "misc") <- NULL
    expect_identical(
        interestingGroups(seurat),
        "sampleName"
    )
    expect_error(
        interestingGroups(seurat) <- "sampleName",
        "seurat object was not generated with bcbioSingleCell"
    )
})

test_that("Assignment method", {
    interestingGroups(bcb) <- "sampleName"
    expect_identical(
        interestingGroups(bcb),
        "sampleName"
    )

    interestingGroups(seurat) <- "sampleName"
    expect_identical(
        interestingGroups(seurat),
        "sampleName"
    )
})
