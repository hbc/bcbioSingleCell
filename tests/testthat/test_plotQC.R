context("plotQC")

test_that("grid", {
    # Example dataset doesn't have a cellular barcode cutoff because we removed
    # the bcbio commands log file (which conflicts with Travis CI)
    p <- plotQC(bcb, return = "grid")
    expect_is(p, "ggplot")
})

test_that("list", {
    p <- plotQC(bcb, return = "list")
    expect_identical(
        names(p),
        c("plotReadsPerCell",
          "plotCellCounts",
          "plotUMIsPerCell",
          "plotGenesPerCell",
          "plotUMIsVsGenes",
          "plotMitoRatio",
          "plotNovelty")
    )
})

test_that("markdown", {
    output <- capture.output(plotQC(bcb, return = "markdown"))
    sep <- c("", "", "")
    expect_identical(
        head(output, 3L),
        c("", "", "## Filtered quality control metrics {.tabset}")
    )
})

test_that("seurat", {
    p <- plotQC(seurat)
    expect_is(p, "ggplot")
})
