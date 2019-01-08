context("Plots")

data(indrops, envir = environment())



# plotQC =======================================================================
test_that("plotQC", {
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

    # Markdown.
    output <- capture.output(plotQC(indrops, return = "markdown"))
    sep <- c("", "", "")
    expect_identical(
        head(output, 3L),
        c("", "", "## Quality control metrics {.tabset}")
    )
})



# plotReadsPerCell =============================================================
# Example dataset doesn't have a cellular barcode cutoff because we removed the
# bcbio commands log file (which conflicts with Travis CI).
test_that("plotReadsPerCell", {
    # Histogram
    x <- plotReadsPerCell(indrops, geom = "histogram")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomStep"
    )

    # Ridgeline
    x <- plotReadsPerCell(indrops, geom = "ridgeline")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomDensityRidges"
    )

    # Violin
    x <- plotReadsPerCell(indrops, geom = "violin")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomViolin"
    )

    # ECDF
    x <- plotReadsPerCell(indrops, geom = "ecdf")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomStep"
    )
})
