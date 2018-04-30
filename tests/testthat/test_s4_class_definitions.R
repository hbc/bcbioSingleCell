context("S4 Class Definitions")



# bcbioSingleCell ==============================================================
test_that("bcbioSingleCell", {
    uploadDir <- system.file("extdata/indrop", package = "bcbioSingleCell")
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



# gene2symbol ==================================================================
colnames <- c("geneID", "geneName")

test_that("gene2symbol : bcbioSingleCell", {
    x <- gene2symbol(indrops_small)
    expect_is(x, "data.frame")
    expect_identical(colnames(x), colnames)
})

test_that("gene2symbol : seurat", {
    x <- gene2symbol(seurat_small)
    expect_is(x, "data.frame")
    expect_identical(colnames(x), colnames)
})




# interestingGroups ============================================================
test_that("interestingGroups : bcbioSingleCell", {
    expect_identical(
        interestingGroups(indrops_small),
        "sampleName"
    )
})

test_that("interestingGroups<- : bcbioSingleCell", {
    error <- "The interesting groups \"XXX\" are not defined"
    expect_error(
        interestingGroups(indrops_small) <- "XXX",
        error
    )
    expect_error(
        interestingGroups(seurat_small) <- "XXX",
        error
    )
})

test_that("interestingGroups : seurat", {
    expect_identical(
        interestingGroups(seurat_small),
        "sampleName"
    )
    expect_identical(
        interestingGroups(Seurat::pbmc_small),
        NULL
    )
})

test_that("interestingGroups<- : seurat", {
    interestingGroups(indrops_small) <- "sampleName"
    expect_identical(
        interestingGroups(indrops_small),
        "sampleName"
    )
    interestingGroups(seurat_small) <- "sampleName"
    expect_identical(
        interestingGroups(seurat_small),
        "sampleName"
    )
    x <- Seurat::pbmc_small
    expect_error(interestingGroups(x) <- "sampleName")
})



# sampleData ===================================================================
clean <- DataFrame(
    "sampleName" = factor("rep_1"),
    "index" = factor("1"),
    "sequence" = factor("TTTTTTTT"),
    "aggregate" = factor("sample"),
    "revcomp" = factor("AAAAAAAA"),
    "interestingGroups" = factor("rep_1"),
    row.names = factor("multiplexed_AAAAAAAA")
)

all <- list(
    "sampleName"  = "factor",
    "fileName"  = "factor",
    "description"  = "factor",
    "index" = "factor",
    "sequence" = "factor",
    "aggregate" = "factor",
    "revcomp" = "factor",
    "interestingGroups" = "factor"
)

test_that("sampleData : bcbioSingleCell", {
    # Clean mode (factor columns only)
    x <- sampleData(indrops_small, clean = TRUE)
    expect_identical(x, clean)

    # Return all columns
    x <- sampleData(indrops_small, clean = FALSE)
    expect_identical(lapply(x, class), all)
})

test_that("sampleData : seurat", {
    # Clean mode (factor columns only)
    x <- sampleData(seurat_small, clean = TRUE)
    expect_identical(x, clean)

    # Return all columns
    x <- sampleData(seurat_small, clean = FALSE)
    expect_identical(lapply(x, class), all)

    # Minimal data for other seurat objects
    x <- sampleData(Seurat::pbmc_small)
    y <- DataFrame(
        "sampleName" = factor("SeuratProject"),
        "interestingGroups" = factor("SeuratProject"),
        row.names = "SeuratProject"
    )
    expect_identical(x, y)
})
