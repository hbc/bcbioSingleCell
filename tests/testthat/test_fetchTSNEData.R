context("fetchTSNEData")

test_that("fetchTSNEData", {
    x <- fetchTSNEData(pbmc_small)
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            "tSNE1" = "numeric",
            "tSNE2" = "numeric",
            "ident" = "factor",
            "nGene" = "integer",
            "nUMI" = "integer",
            "origIdent" = "factor",
            "res0x8" = "factor",
            "res1" = "factor",
            "centerX" = "numeric",
            "centerY" = "numeric"
        )
    )

    subset <- head(x, 1L) %>%
        rownames_to_column() %>%
        mutate_if(is.numeric, funs(round(., digits = 3L))) %>%
        column_to_rownames()
    target <- data.frame(
        "tSNE1" = 14.422,
        "tSNE2" = 8.336,
        "ident" = factor("0", levels = c("0", "1", "2", "3")),
        "nGene" = 47L,
        "nUMI" = 70L,
        "origIdent" = factor("SeuratProject"),
        "res0x8" = factor("0"),
        "res1" = factor("0", levels = c("0", "1", "2", "3")),
        "centerX" = -0.332,
        "centerY" = 18.76,
        row.names = "ATGCCAGAACGACT",
        stringsAsFactors = FALSE
    )
    expect_equal(subset, target)
})
