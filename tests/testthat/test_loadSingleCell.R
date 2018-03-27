context("loadSingleCell")

test_that("Homo sapiens", {
    uploadDir <- system.file("extdata/indrop", package = "bcbioSingleCell")
    sampleMetadataFile <- file.path(uploadDir, "metadata.csv")
    x <- loadSingleCell(
        uploadDir = uploadDir,
        sampleMetadataFile = sampleMetadataFile,
        organism = "Homo sapiens"
    )
    expect_s4_class(x, "bcbioSingleCell")
})
