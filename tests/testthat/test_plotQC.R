context("plotQC")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("grid", {
    # Example dataset doesn't have a cellular barcode cutoff because we removed
    # the bcbio commands log file (which conflicts with Travis CI)
    p <- suppressWarnings(plotQC(bcb, return = "grid"))
    expect_is(p, "ggplot")
})

test_that("list", {
    p <- suppressWarnings(plotQC(bcb, return = "list"))
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
    output <- capture.output(suppressWarnings(
        plotQC(bcb, return = "markdown")
    ))
    sep <- c("", "", "")
    expect_identical(
        output,
        c("",
          "",
          "## Filtered quality control metrics {.tabset}",
          sep,
          "### Reads per cell",
          sep,
          "### Cell counts",
          sep,
          "### UMI counts per cell",
          sep,
          "### Genes detected",
          sep,
          "### UMIs vs. genes",
          sep,
          "### Mitochondrial counts ratio",
          sep,
          "### Novelty",
          "")
    )
})

test_that("seurat", {
    p <- plotQC(seurat)
    expect_is(p, "ggplot")
})
