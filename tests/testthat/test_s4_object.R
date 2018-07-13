context("S4 Object : bcbioSingleCell")



# bcbioSingleCell ==============================================================
test_that("bcbioSingleCell", {
    uploadDir <- system.file("extdata/indrops", package = "bcbioSingleCell")
    sampleMetadataFile <- file.path(uploadDir, "metadata.csv")

    # Organism
    x <- bcbioSingleCell(
        uploadDir = uploadDir,
        sampleMetadataFile = sampleMetadataFile,
        organism = "Homo sapiens"
    )
    expect_s4_class(x, "bcbioSingleCell")

    # NULL organism
    x <- suppressWarnings(
        bcbioSingleCell(
            uploadDir = uploadDir,
            sampleMetadataFile = sampleMetadataFile,
            organism = NULL
        )
    )
    expect_s4_class(x, "bcbioSingleCell")
})
