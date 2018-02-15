context("interestingGroups")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("bcbioSingleCell", {
    expect_identical(
        interestingGroups(bcb),
        "sampleName"
    )
})

test_that("Assignment of undefined interesting group", {
    error <- paste(
        "is_subset : The element 'XXX' in interestingGroups is not",
        "in colnames\\(object\\)"
    )
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
    expect_identical(
        interestingGroups(seurat),
        "sampleName"
    )
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
