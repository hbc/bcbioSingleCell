context("plotQC")

test_that("grid", {
    object <- indrops

    # Grid. Example dataset doesn't have a cellular barcode cutoff because we
    # removed the bcbio commands log file (which conflicts with Travis CI).
    p <- plotQC(object, return = "grid")
    expect_s3_class(p, "ggplot")

    # List. "Cell Counts" only returns for filtered datasets. "Reads per Cell"
    # only returns for bcbioSingleCell.
    x <- plotQC(object, return = "list")
    names <- c(
        "genesPerCell",
        "mitoRatio",
        "novelty",
        "umisPerCell",
        "umisVsGenes",
        "zerosVsDepth"
    )
    expect_true(all(names %in% names(x)))
})

test_that("markdown", {
    output <- capture.output(plotQC(indrops, return = "markdown"))
    sep <- c("", "", "")
    expect_identical(
        head(output, 3L),
        c("", "", "## Quality control metrics {.tabset}")
    )
})
