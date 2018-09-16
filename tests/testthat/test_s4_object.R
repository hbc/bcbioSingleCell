context("S4 Objects")



# bcbioSingleCell ==============================================================
# Need to add YAML metadata to test without `sampleMetadataFile`
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
    x <- bcbioSingleCell(
        uploadDir = uploadDir,
        sampleMetadataFile = sampleMetadataFile,
        organism = NULL
    )
    expect_s4_class(x, "bcbioSingleCell")
})



# CellRanger ===================================================================
test_that("CellRanger", {
    object <- suppressWarnings(CellRanger(
        uploadDir = system.file(
            "extdata/cellranger",
            package = "bcbioSingleCell"
        ),
        organism = "Homo sapiens"
    ))
    expect_is(object, "CellRanger")
    expect_identical(dim(object), c(500L, 500L))
    expect_identical(
        sampleNames(object),
        c(pbmc_1 = "pbmc")
    )
    expect_identical(
        object = head(colnames(object), n = 1L),
        expected = "pbmc_1_AAACCTGAGAAGGCCT"
    )
    expect_identical(
        object = head(rownames(object), n = 1L),
        expected = "ENSG00000004487"
    )
})
