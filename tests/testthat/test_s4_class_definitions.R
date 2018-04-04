context("S4 Class Definitions")



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
        interestingGroups(pbmc_small),
        "sampleName"
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
    expect_error(
        interestingGroups(pbmc_small) <- "sampleName",
        "object was not created with bcbioSingleCell"
    )
})



# sampleData ===================================================================
target <- DataFrame(
    "sampleID" = factor("multiplexed_AAAAAAAA"),
    "sampleName" = factor("rep_1"),
    "description" = factor("multiplexed"),
    "fileName" = factor("multiplexed.fastq.gz"),
    "index" = factor("1"),
    "sequence" = factor("TTTTTTTT"),
    "sampleNameAggregate" = factor("sample"),
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
    expect_identical(x, target)
})

test_that("sampleData : seurat", {
    x <- sampleData(seurat_small)
    expect_identical(x, target)

    x <- sampleData(pbmc_small)
    y <- DataFrame(
        "sampleID" = factor("SeuratProject"),
        "sampleName" = factor("SeuratProject"),
        "description" = factor("SeuratProject"),
        "interestingGroups" = factor("SeuratProject"),
        row.names = "SeuratProject"
    )
    expect_identical(x, y)
})
