context("interestingGroups")

test_that("bcbioSingleCell", {
    expect_identical(
        interestingGroups(bcb_small),
        "sampleName"
    )
})

test_that("Assignment of undefined interesting group", {
    error <- paste(
        "is_subset : The element 'XXX' in interestingGroups is not",
        "in colnames\\(x\\)"
    )
    expect_error(
        interestingGroups(bcb_small) <- "XXX",
        error
    )
    expect_error(
        interestingGroups(seurat_small) <- "XXX",
        error
    )
})

test_that("seurat", {
    expect_identical(
        interestingGroups(seurat_small),
        "sampleName"
    )
    expect_identical(
        interestingGroups(pbmc_small),
        "sampleName"
    )
})

test_that("Assignment method", {
    interestingGroups(bcb_small) <- "sampleName"
    expect_identical(
        interestingGroups(bcb_small),
        "sampleName"
    )
    interestingGroups(seurat_small) <- "sampleName"
    expect_identical(
        interestingGroups(seurat_small),
        "sampleName"
    )
    expect_error(
        interestingGroups(pbmc_small) <- "sampleName",
        "object was not created with bcbioSingleCell"
    )
})
