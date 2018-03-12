context("loadSingleCell")

test_that("Homo sapiens", {
    extdataDir <- system.file("extdata", package = "bcbioSingleCell")
    uploadDir <- file.path(extdataDir, "harvard_indrop_v3")
    sampleMetadataFile <- file.path(extdataDir, "harvard_indrop_v3.xlsx")
    bcb <- suppressWarnings(suppressMessages(
        loadSingleCell(
        uploadDir = uploadDir,
        sampleMetadataFile = sampleMetadataFile,
        organism = "Homo sapiens")
    ))
    expect_is(bcb, "bcbioSingleCell")
})
