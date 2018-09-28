context("S4 Objects")



# bcbioSingleCell ==============================================================
test_that("bcbioSingleCell", {
    uploadDir <- system.file("extdata/indrops", package = "bcbioSingleCell")

    # Minimal mode, with no metadata or annotations.
    # This is fast but doesn't slot a lot of useful info.
    x <- bcbioSingleCell(uploadDir = uploadDir)
    expect_s4_class(x, "bcbioSingleCell")

    # User-defined metadata.
    x <- bcbioSingleCell(
        uploadDir = uploadDir,
        sampleMetadataFile <- file.path(uploadDir, "metadata.csv")
    )
    expect_s4_class(x, "bcbioSingleCell")

    # Automatic organism annotations from AnnotationHub.
    x <- bcbioSingleCell(
        uploadDir = uploadDir,
        organism = "Homo sapiens"
    )
    expect_s4_class(x, "bcbioSingleCell")
})



# CellRanger ===================================================================
test_that("CellRanger", {
    uploadDir = system.file("extdata/cellranger", package = "bcbioSingleCell")

    # Minimal mode, with no metadata or annotations.
    # This is fast but doesn't slot a lot of useful info.
    x <- CellRanger(uploadDir = uploadDir)
    expect_s4_class(x, "CellRanger")

    # Automatic organism annotations from AnnotationHub.
    x <- CellRanger(
        uploadDir = uploadDir,
        organism = "Homo sapiens",
        ensemblRelease = 87L
    )
    expect_s4_class(x, "CellRanger")

    # Minimal example contains some genes that are dead on current Ensembl.
    expect_warning(
        object = CellRanger(
            uploadDir = uploadDir,
            organism = "Homo sapiens",
            ensemblRelease = 92L
        ),
        regexp = "Genes missing in rowRanges."
    )
})
