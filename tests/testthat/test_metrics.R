context("metrics")

test_that("metrics : bcbioSingleCell", {
    x <- metrics(bcb_small)
    expect_identical(
        lapply(x, class),
        list(
            "nCount" = "integer",
            "nUMI" = "integer",
            "nGene" = "integer",
            "nCoding" = "integer",
            "nMito" = "integer",
            "log10GenesPerUMI" = "numeric",
            "mitoRatio" = "numeric",
            "sampleID" = "factor",
            "sampleName" = "factor",
            "description" = "factor",
            "fileName" = "factor",
            "index" = "factor",
            "sequence" = "factor",
            "sampleNameAggregate" = "factor",
            "revcomp" = "factor",
            "interestingGroups" = "factor"
        )
    )
})

test_that("metrics : seurat", {
    x <- metrics(pbmc_small)
    expect_identical(
        lapply(x, class),
        list(
            "nGene" = "integer",
            "nUMI" = "integer",
            "origIdent" = "factor",
            "res0x8" = "factor",
            "res1" = "factor",
            "sampleID" = "factor",
            "sampleName" = "factor",
            "description" = "factor",
            "interestingGroups" = "factor",
            "ident" = "factor"
        )
    )
})
