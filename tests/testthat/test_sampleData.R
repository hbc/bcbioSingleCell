context("sampleData")

target <- data.frame(
    "sampleID" = "multiplexed_AAAAAAAA",
    "sampleName" = "rep_1",
    "description" = "multiplexed",
    "fileName" = "multiplexed.fastq.gz",
    "index" = "1",
    "sequence" = "TTTTTTTT",
    "sampleNameAggregate" = "sample",
    "revcomp" = "AAAAAAAA",
    "interestingGroups" = "rep_1",
    row.names = "multiplexed_AAAAAAAA",
    stringsAsFactors = TRUE
)

test_that("sampleData : bcbioSingleCell", {
    x <- sampleData(bcb_small)
    expect_is(x, "data.frame")
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
    y <- data.frame(
        "sampleID" = "SeuratProject",
        "sampleName" = "SeuratProject",
        "description" = "SeuratProject",
        "interestingGroups" = "SeuratProject",
        row.names = "SeuratProject",
        stringsAsFactors = TRUE
    )
    expect_identical(x, y)
})
