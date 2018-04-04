context("sampleData")

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
