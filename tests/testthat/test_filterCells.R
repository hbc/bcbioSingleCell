context("filterCells")

test_that("Default parameters", {
    filtered <- filterCells(bcb)
    expect_is(filtered, "bcbioSingleCell")
    expect_identical(
        dim(filtered),
        c(1000L, 421L)
    )
    expect_is(
        metadata(filtered)[["filterParams"]],
        "list"
    )
    expect_is(
        metadata(filtered)[["filterCells"]],
        "character"
    )
    expect_is(
        metadata(filtered)[["filterGenes"]],
        "character"
    )
    expect_identical(
        metadata(filtered)[["subset"]],
        TRUE
    )
})

test_that("Capture output", {
    output <- capture.output(filterCells(bcb))
    expect_is(output, "character")
    expect_identical(
        output[[1L]],
        "Filtering parameters:"
    )
})

test_that("Maximum parameters", {
    # This should return an object with the same dimensions
    filtered <- filterCells(
        bcb,
        minUMIs = 0L,
        maxUMIs = Inf,
        minGenes = 0L,
        maxGenes = Inf,
        maxMitoRatio = 1L,
        minNovelty = 0L,
        minCellsPerGene = 0L)
    expect_is(filtered, "bcbioSingleCell")
    expect_identical(
        dim(filtered),
        dim(bcb)
    )
})

test_that("Cutoff failures", {
    expect_error(
        filterCells(bcb, minUMIs = Inf),
        "No cells passed `minUMIs` cutoff"
    )
    expect_error(
        filterCells(bcb, minCellsPerGene = Inf),
        "No genes passed `minCellsPerGene` cutoff"
    )
})

test_that("Per sample cutoffs", {
    # Get the count of sample1 (run1_AGAGGATA)
    # We're applying no filtering to that sample
    metrics <- metrics(bcb)
    ncells1 <- length(which(metrics[["sampleID"]] == "run1_AGAGGATA"))
    filtered <- filterCells(
        bcb,
        minUMIs = c(
            "run1_AGAGGATA" = 0L,
            "run2_AGAGGATA" = Inf),
        maxUMIs = c(
            "run1_AGAGGATA" = Inf,
            "run2_AGAGGATA" = 0L),
        minGenes = c(
            "run1_AGAGGATA" = 0L,
            "run2_AGAGGATA" = Inf),
        maxGenes = c(
            "run1_AGAGGATA" = Inf,
            "run2_AGAGGATA" = 0L),
        maxMitoRatio = c(
            "run1_AGAGGATA" = 0L,
            "run2_AGAGGATA" = Inf),
        minNovelty = c(
            "run1_AGAGGATA" = 0L,
            "run2_AGAGGATA" = Inf)
    )
    expect_identical(
        ncol(filtered),
        ncells1
    )
})
