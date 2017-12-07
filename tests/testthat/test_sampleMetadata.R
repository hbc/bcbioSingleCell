context("sampleMetadata")

test_that("sampleMetadata", {
    data <- sampleMetadata(bcb)
    expect_is(data, "data.frame")
    expect_identical(
        lapply(data, class),
        list(
            "sampleID" = "factor",
            "sampleName"  = "factor",
            "description" = "factor",
            "fileName" = "factor",
            "index" = "factor",
            "sequence" = "factor",
            "sampleNameAggregate" = "factor",
            "sequencingReplicate" = "factor",
            "biologicalReplicate" = "factor",
            "genotype" = "factor",
            "concNgUl" = "factor",
            "revcomp" = "factor",
            "interestingGroups" = "factor"
        )
    )
})
