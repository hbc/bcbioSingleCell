context("Plot Functions")



# plotQC =======================================================================
with_parameters_test_that(
    "plotQC", {
        # Grid
        # Example dataset doesn't have a cellular barcode cutoff because we
        # removed the bcbio commands log file (which conflicts with Travis CI).
        p <- plotQC(object, return = "grid")
        expect_is(p, "ggplot")

        # List
        x <- plotQC(object, return = "list")
        # "Cell Counts" only returns for filtered datasets.
        # "Reads per Cell" only returns for bcbioSingleCell.
        names <- c(
            "genesPerCell",
            "mitoRatio",
            "novelty",
            "umisPerCell",
            "umisVsGenes",
            "zerosVsDepth"
        )
        expect_true(all(names %in% names(x)))

        # Markdown
        output <- capture.output(plotQC(indrops_small, return = "markdown"))
        sep <- c("", "", "")
        expect_identical(
            head(output, 3L),
            c("", "", "## Quality control metrics {.tabset}")
        )
    },
    object = list(
        bcbioSingleCell = indrops_small,
        CellRanger = cellranger_small
    )
)



# plotReadsPerCell =============================================================
# Example dataset doesn't have a cellular barcode cutoff because we removed the
# bcbio commands log file (which conflicts with Travis CI).
test_that("plotReadsPerCell", {
    # Histogram
    x <- plotReadsPerCell(indrops_small, geom = "histogram")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomStep"
    )

    # Ridgeline
    x <- plotReadsPerCell(indrops_small, geom = "ridgeline")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomDensityRidges"
    )

    # Violin
    x <- plotReadsPerCell(indrops_small, geom = "violin")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomViolin"
    )

    # ECDF
    x <- plotReadsPerCell(indrops_small, geom = "ecdf")
    expect_is(x, "ggplot")
    expect_is(
        x %>%
            .[["layers"]] %>%
            .[[1L]] %>%
            .[["geom"]],
        "GeomStep"
    )
})
