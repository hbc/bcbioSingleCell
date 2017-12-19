context("fetchTSNEExpressionData")

load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

symbol <- counts(seurat) %>% rownames() %>% head()
ensgene <- bcbio(seurat, "gene2symbol") %>%
    .[which(.[["symbol"]] %in% symbol), "ensgene", drop = TRUE]

test_that("symbol", {
    data <- fetchTSNEExpressionData(seurat, genes = symbol, format = "symbol")
    expect_is(data, "grouped_df")
    group <- group_vars(data)
    expect_identical(group, "gene")
    expect_identical(
        lapply(data, class),
        list(
            "gene" = "character",
            "cellID" = "character",
            "expression" = "numeric",
            "geomean" = "numeric",
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
        mutate_if(is.numeric, funs(round(., digits = 3L)))
    target <- tibble(
        "gene" = "CD99",
        "cellID" = "M1_AAACACTA_CTTAGGTA",
        "expression" = 1.762,
        "geomean" = 1.165,
        "tSNE1" = -3.056,
        "tSNE2" = -4.304,
        "ident" = factor("0", levels = c("0", "1", "2", "3")),
        "nGene" = 734,  # numeric here
        "nUMI" = 6220,
        "sampleID" = factor("M1"),
        "nCoding" = 5495,
        "nMito" = 442,
        "log10GenesPerUMI" = 0.755,
        "mitoRatio" = 0.071,
        "nCount" = 95155,  # numeric here
        "sampleName" = factor("M1"),
        "description" = factor("M1"),
        "interestingGroups" = factor("M1"),
        "orig.ident" = factor("M1"),
        "res.0.8" = "0",
        "centerX" = -6.651,  # reorder
        "centerY" = -4.226  # reorder
    ) %>%
        group_by(!!sym(group))
    expect_equal(subset, target)
})
