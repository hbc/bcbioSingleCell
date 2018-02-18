context("fetchTSNEExpressionData")

load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))

genes <- counts(seurat) %>% rownames() %>% head()

test_that("seurat", {
    data <- fetchTSNEExpressionData(seurat, genes = genes)
    expect_is(data, "data.frame")
    expect_identical(
        lapply(data, class),
        list(
            "tSNE1" = "numeric",
            "tSNE2" = "numeric",
            "ident" = "factor",
            "nGene" = "integer",
            "nUMI" = "integer",
            "sampleID" = "factor",
            "nCoding" = "integer",
            "nMito" = "integer",
            "log10GenesPerUMI" = "numeric",
            "mitoRatio" = "numeric",
            "nCount" = "integer",
            "sampleName" = "factor",
            "description" = "factor",
            "interestingGroups" = "factor",
            "origIdent" = "factor",
            "res0x8" = "factor",
            "centerX" = "numeric",
            "centerY" = "numeric",
            "mean" = "numeric",
            "median" = "numeric",
            "sum" = "numeric"  # FIXME coerce to integer
        )
    )

    subset <- data[1L, ] %>%
        rownames_to_column() %>%
        mutate_if(is.numeric, funs(round(., digits = 3L))) %>%
        column_to_rownames()
    # The round function will coerce integers to numerics. This is the rationale
    # for the `as.numeric()` usage below.
    identLevels <- c("0", "1", "2", "3")
    target <- data.frame(
        "tSNE1" = -3.056,
        "tSNE2" = -4.304,
        "ident" = factor("0", levels = identLevels),
        "nGene" = as.numeric(734L),
        "nUMI" = as.numeric(6220L),
        "sampleID" = factor("M1"),
        "nCoding" = as.numeric(5495L),
        "nMito" = as.numeric(442L),
        "log10GenesPerUMI" = 0.755,
        "mitoRatio" = 0.071,
        "nCount" = as.numeric(95155L),
        "sampleName" = factor("M1"),
        "description" = factor("M1"),
        "interestingGroups" = factor("M1"),
        "origIdent" = factor("M1"),
        "res0x8" = factor("0", levels = identLevels),
        "centerX" = -6.651,
        "centerY" = -4.226,
        "mean" = 1.165,
        "median" = 1.36,
        "sum" = 6.988,
        row.names = "M1_AAACACTA_CTTAGGTA",
        stringsAsFactors = TRUE
    )
    expect_identical(subset, target)
})
