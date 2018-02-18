context("sampleMetadata")

load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))

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

seuratdata <- data.frame(
    "sampleID" = "M1",
    "sampleName" = "M1",
    "description" = "M1",
    "interestingGroups" = "M1",
    row.names = "M1",
    stringsAsFactors = TRUE)

test_that("seurat", {
    data <- sampleMetadata(seurat)
    expect_is(data, "data.frame")
    expect_identical(data, seuratdata)
})

test_that("seurat without stashed metadata", {
    bcbio(seurat, "sampleMetadata") <- NULL
    data <- sampleMetadata(seurat)
    expect_is(data, "data.frame")
    expect_identical(data, seuratdata)
})
