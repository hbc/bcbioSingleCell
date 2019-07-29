context("calculateMetrics")

## FIXME counts
## FIXME rowRanges
## FIXME prefilter

## FIXME calculateMetrics rowRanges setdiff
## FIXME calculateMetrics with broadClass
## FIXME calculateMetrics without broadClass

counts <- assay(indrops)

test_that("No biotype information from rowRanges", {
    x <- calculateMetrics(counts = counts)
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
})
