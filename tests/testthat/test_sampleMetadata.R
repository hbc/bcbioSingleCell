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

test_that("prepareSampleMetadataFromSeurat", {
    meta <- slot(seurat, "meta.data")
    expect_is(meta, "data.frame")
    data <- .prepareSampleMetadataFromSeurat(meta)
    expect_is(data, "data.frame")
    expect_identical(
        data,
        data.frame(
            "sampleID" = "M1",
            "sampleName" = "M1",
            "description" = "M1",
            "interestingGroups" = "M1",
            "origIdent" = "M1",
            stringsAsFactors = TRUE)
    )
})
