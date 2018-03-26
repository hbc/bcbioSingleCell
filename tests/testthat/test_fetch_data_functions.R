context("Fetch Data Functions")



# fetchPCAData =================================================================
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
    target <- data.frame(
        "pc1" = 0.443,
        "pc2" = -1.575,
        "ident" = factor("0", levels = c("0", "1", "2", "3")),
        "nGene" = 47L,
        "nUMI" = 70L,
        "origIdent" = factor("SeuratProject"),
        "res0x8" = factor("0", levels = "0"),
        "res1" = factor("0", levels = c("0", "1", "2", "3")),
        "centerX" = -1.176,
        "centerY" = -1.683,
        row.names = "ATGCCAGAACGACT",
        stringsAsFactors = FALSE
    )
    # Use `glimpse()` to compare the rounded numbers
    expect_equal(subset, target)
})



# fetchTSNEData ================================================================
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
        tibble::rownames_to_column(.) %>%
        dplyr::mutate_if(
            is.numeric,
            dplyr::funs(round(., digits = 3L))
        ) %>%
        tibble::column_to_rownames(.)
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



# fetchTSNEExpressionData ======================================================
test_that("fetchTSNEExpressionData", {
    x <- fetchTSNEExpressionData(
        object = pbmc_small,
        genes = head(rownames(pbmc_small))
    )
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
            "centerY" = "numeric",
            "mean" = "numeric",
            "median" = "numeric",
            "sum" = "numeric"
        )
    )

    subset <- x[1L, ] %>%
        tibble::rownames_to_column(.) %>%
        dplyr::mutate_if(is.numeric, dplyr::funs(round(., digits = 3L))) %>%
        tibble::column_to_rownames(.)
    # The round function will coerce integers to numerics. This is the rationale
    # for the `as.numeric()` usage below.
    target <- data.frame(
        "tSNE1" = 14.422,
        "tSNE2" = 8.336,
        "ident" = factor("0", levels = c("0", "1", "2", "3")),
        "nGene" = as.numeric(47L),
        "nUMI" = as.numeric(70L),
        "origIdent" = factor("SeuratProject"),
        "res0x8" = factor("0"),
        "res1" = factor("0", levels = c("0", "1", "2", "3")),
        "centerX" = -0.332,
        "centerY" = 18.76,
        "mean" = 1.656,
        "median" = 0L,
        "sum" = 9.938,
        row.names = "ATGCCAGAACGACT",
        stringsAsFactors = TRUE
    )
    expect_equal(subset, target)
})
