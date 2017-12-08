context("sampleMetadata")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("bcbioSingleCell", {
    data <- sampleMetadata(bcb)
    expect_is(data, "data.frame")
    expect_identical(
        lapply(data, class),
        list(
            "sampleID" = "factor",
            "sampleName"  = "factor",
            "description" = "factor",
            "fileName" = "factor",
            "index" = "factor",
            "sequence" = "factor",
            "sampleNameAggregate" = "factor",
            "revcomp" = "factor",
            "interestingGroups" = "factor"
        )
    )
})

test_that("seurat", {
    data <- sampleMetadata(seurat)
    expect_is(data, "data.frame")
    expect_identical(
        lapply(data, class),
        list(
            "sampleID" = "factor",
            "sampleName"  = "factor",
            "description" = "factor",
            "interestingGroups" = "factor"
        )
    )
})
