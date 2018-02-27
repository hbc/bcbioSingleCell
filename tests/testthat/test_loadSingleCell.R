context("loadSingleCell")

test_that("Homo sapiens", {
    extdataDir <- system.file("extdata", package = "bcbioSingleCell")
    uploadDir <- file.path(extdataDir, "harvard_indrop_v3")
    sampleMetadataFile <- file.path(extdataDir, "harvard_indrop_v3.xlsx")
    bcb <- loadSingleCell(
        uploadDir = uploadDir,
        sampleMetadataFile = sampleMetadataFile)
    expect_is(bcb, "bcbioSingleCell")
})
