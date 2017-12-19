context("fetchTSNEData")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("fetchTSNEData", {
    data <- fetchTSNEData(seurat)
    expect_is(data, "data.frame")
    expect_identical(
        lapply(data, class),
        list(
            "tSNE1" = "numeric",
            "tSNE2" = "numeric",
            "ident" = "factor",
            "nGene" = "integer",
            "nUMI" = "numeric",  # FIXME integer
            "sampleID" = "factor",
            "nCoding" = "numeric",  # FIXME integer
            "nMito" = "numeric",  # FIXME integer
            "log10GenesPerUMI" = "numeric",
            "mitoRatio" = "numeric",
            "nCount" = "integer",
            "sampleName" = "factor",
            "description" = "factor",
            "interestingGroups" = "factor",
            "orig.ident" = "factor",  # camel?
            "res.0.8" = "character",  # camel?
            "centerX" = "numeric",
            "centerY" = "numeric"
        )
    )
    subset <- data[1L, ] %>%
        rownames_to_column() %>%
        mutate_if(is.numeric, funs(round(., digits = 3L))) %>%
        column_to_rownames()
    target <- data.frame(
        "tSNE1" = -3.056,
        "tSNE2" = -4.304,
        "ident" = factor("0", levels = c("0", "1", "2", "3")),
        "nGene" = 734L,
        "nUMI" = 6220,
        "sampleID" = factor("M1"),
        "nCoding" = 5495,
        "nMito" = 442,
        "log10GenesPerUMI" = 0.755,
        "mitoRatio" = 0.071,
        "nCount" = 95155L,
        "sampleName" = factor("M1"),
        "description" = factor("M1"),
        "interestingGroups" = factor("M1"),
        "orig.ident" = factor("M1"),
        "res.0.8" = "0",
        "centerX" = -6.651,
        "centerY" = -4.226,
        row.names = "M1_AAACACTA_CTTAGGTA",
        stringsAsFactors = FALSE
    )
    expect_equal(subset, target)
})
