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
    x <- gene2symbol(bcb_small)
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
        interestingGroups(bcb_small),
        "sampleName"
    )
})

test_that("interestingGroups<- : bcbioSingleCell", {
    error <- paste(
        "is_subset : The element 'XXX' in interestingGroups is not",
        "in colnames\\(x\\)"
    )
    expect_error(
        interestingGroups(bcb_small) <- "XXX",
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
    interestingGroups(bcb_small) <- "sampleName"
    expect_identical(
        interestingGroups(bcb_small),
        "sampleName"
    )
    interestingGroups(seurat_small) <- "sampleName"
    expect_identical(
        interestingGroups(seurat_small),
        "sampleName"
    )
    x <- Seurat::pbmc_small
    expect_error(
        interestingGroups(x) <- "sampleName",
        "object was not created with bcbioSingleCell"
    )
})



# sampleData ===================================================================
target <- DataFrame(
    "sampleName" = factor("rep_1"),
    "fileName" = factor("multiplexed.fastq.gz"),
    "description" = factor("multiplexed"),
    "index" = factor("1"),
    "sequence" = factor("TTTTTTTT"),
    "aggregate" = factor("sample"),
    "revcomp" = factor("AAAAAAAA"),
    "interestingGroups" = factor("rep_1"),
    row.names = factor("multiplexed_AAAAAAAA")
)

test_that("sampleData : bcbioSingleCell", {
    x <- sampleData(bcb_small)
    expect_is(x, "DataFrame")
    expect_identical(
        lapply(x, class),
        list(
            "sampleName"  = "factor",
            "fileName" = "factor",
            "description" = "factor",
            "index" = "factor",
            "sequence" = "factor",
            "aggregate" = "factor",
            "revcomp" = "factor",
            "interestingGroups" = "factor"
        )
    )
    expect_identical(x, target)
})

test_that("sampleData : seurat", {
    x <- sampleData(seurat_small)
    expect_identical(x, target)

    x <- sampleData(Seurat::pbmc_small)
    y <- DataFrame(
        "sampleName" = factor("SeuratProject"),
        "interestingGroups" = factor("SeuratProject"),
        row.names = "SeuratProject"
    )
    expect_identical(x, y)
})
