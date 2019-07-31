context("calculateMetrics")

counts <- counts(indrops)
rowRanges <- rowRanges(indrops)

test_that("rowRanges argument", {
    x <- calculateMetrics(counts = counts, rowRanges = NULL)
    expect_identical(
        object = colnames(x),
        expected = c(
            "nUMI",
            "nGene",
            "nCoding",
            "nMito",
            "log10GenesPerUMI",
            "mitoRatio"
        )
    )
    expect_identical(
        object = vapply(X = x, FUN = anyNA, FUN.VALUE = logical(1L)),
        expected = c(
            nUMI = FALSE,
            nGene = FALSE,
            nCoding = TRUE,
            nMito = TRUE,
            log10GenesPerUMI = FALSE,
            mitoRatio = TRUE
        )
    )
    x <- calculateMetrics(counts = counts, rowRanges = rowRanges)
    expect_false(any(vapply(X = x, FUN = anyNA, FUN.VALUE = logical(1L))))
})

test_that("Low pass prefiltering", {
    ## All barcodes in example should pass.
    x <- calculateMetrics(counts = counts, prefilter = TRUE)
    expect_identical(nrow(x), ncol(counts))
    ## Simulate a poor barcode.
    bad <- counts
    bad[, seq_len(2L)] <- 0L
    x <- calculateMetrics(counts = bad, prefilter = TRUE)
    expect_identical(nrow(x), ncol(counts) - 2L)
})
