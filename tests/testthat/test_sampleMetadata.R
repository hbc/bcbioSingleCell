context("sampleMetadata")

test_that("bcbioSingleCell", {
    data <- sampleMetadata(bcb_small)
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

test_that("seurat_small", {
    data <- sampleMetadata(seurat_small)
    expect_is(data, "data.frame")
    expect_identical(data, seuratdata)
})

test_that("seurat_small without stashed metadata", {
    bcbio(seurat_small, "sampleMetadata") <- NULL
    data <- sampleMetadata(seurat_small)
    expect_is(data, "data.frame")
    expect_identical(data, seuratdata)
})
