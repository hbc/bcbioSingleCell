context("Plot Functions")



# plotQC =======================================================================
test_that("plotQC : grid", {
    # Example dataset doesn't have a cellular barcode cutoff because we removed
    # the bcbio commands log file (which conflicts with Travis CI)
    p <- plotQC(indrops_small, return = "grid")
    expect_is(p, "ggplot")
})

test_that("plotQC : list", {
    p <- plotQC(indrops_small, return = "list")
    expect_identical(
        names(p),
        c(
            "Cell Counts",
            "Reads per Cell",
            "UMIs per Cell",
            "Genes per Cell",
            "UMIs vs. Genes",
            "Novelty",
            "Mito Ratio",
            "Zeros vs. Depth"
        )
    )
})

test_that("plotQC : markdown", {
    output <- capture.output(plotQC(indrops_small, return = "markdown"))
    sep <- c("", "", "")
    expect_identical(
        head(output, 3L),
        c("", "", "## Filtered quality control metrics {.tabset}")
    )
})



# plotReadsPerCell =============================================================
# Example dataset doesn't have a cellular barcode cutoff because we removed the
# bcbio commands log file (which conflicts with Travis CI)
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
