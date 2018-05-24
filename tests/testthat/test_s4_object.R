# FIXME Add seurat coercion check
# FIXME Add `updateObject()` test



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



context("S4 Object : seurat")



# bcbio ========================================================================
test_that("bcbio", {
    x <- bcbio(seurat_small)
    expect_is(x, "list")
    expect_identical(
        lapply(x, class),
        list(
            rowRanges = structure("GRanges", package = "GenomicRanges"),
            metadata = "list"
        )
    )
})

test_that("bcbio assignment", {
    x <- seurat_small
    # Stash as new slot
    bcbio(x, "stash") <- "XXX"
    expect_identical(
        bcbio(x, "stash"),
        "XXX"
    )
    # Metadata stash
    bcbio(x, "metadata")[["stash"]] <- "YYY"
    expect_identical(
        bcbio(x, "metadata")[["stash"]],
        "YYY"
    )
})
