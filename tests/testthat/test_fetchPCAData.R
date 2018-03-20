context("fetchPCAData")

test_that("fetchPCAData", {
    x <- fetchPCAData(pbmc_small)
    expect_is(x, "data.frame")
    expect_identical(
        lapply(x, class),
        list(
            "pc1" = "numeric",
            "pc2" = "numeric",
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
    subset <- x[1L, ] %>%
        tibble::rownames_to_column(.) %>%
        dplyr::mutate_if(
            is.numeric,
            dplyr::funs(round(., digits = 3L))
        ) %>%
        tibble::column_to_rownames(.)
    levels <- c("0", "1", "2", "3")
    target <- data.frame(
        "pc1" = 0.443,
        "pc2" = -1.575,
        "ident" = factor("0", levels = levels),
        "nGene" = 47L,
        "nUMI" = 70L,
        "origIdent" = factor("SeuratProject"),
        "res0x8" = factor("0", levels = "0"),
        "res1" = factor("0", levels = levels),
        "centerX" = -1.176,
        "centerY" = -1.683,
        row.names = "ATGCCAGAACGACT",
        stringsAsFactors = FALSE
    )
    # Use `glimpse()` to compare the rounded numbers
    expect_equal(subset, target)
})
