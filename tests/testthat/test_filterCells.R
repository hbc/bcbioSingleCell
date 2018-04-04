context("filterCells")

test_that("filterCells : Default parameters", {
    x <- filterCells(bcb_small, minGenes = 0L)
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(
        dim(x),
        c(500L, 132L)
    )
    expect_is(
        metadata(x)[["filterParams"]],
        "list"
    )
    expect_is(
        metadata(x)[["filterCells"]],
        "character"
    )
    expect_is(
        metadata(x)[["filterGenes"]],
        "character"
    )
    expect_identical(
        metadata(x)[["subset"]],
        TRUE
    )
})

test_that("filterCells : Maximum parameters", {
    # This should return an object with the same dimensions
    x <- filterCells(
        bcb_small,
        minUMIs = 0L,
        maxUMIs = Inf,
        minGenes = 0L,
        maxGenes = Inf,
        maxMitoRatio = 1L,
        minNovelty = 0L,
        minCellsPerGene = 0L
    )
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(
        dim(x),
        dim(bcb_small)
    )
})

test_that("filterCells : Cutoff failures", {
    expect_error(
        filterCells(bcb_small, minUMIs = Inf),
        "No cells passed `minUMIs` cutoff"
    )
})

test_that("filterCells : Per sample cutoffs", {
    # Get the count of sample1 (run1_AGAGGATA)
    # We're applying no filtering to that sample
    metrics <- metrics(bcb_small)
    sample <- levels(metrics[["sampleID"]])[[1L]]
    expect_identical(sample, "multiplexed_AAAAAAAA")
    nCells <- length(which(metrics[["sampleID"]] == sample))
    x <- filterCells(
        object = bcb_small,
        minUMIs = c("multiplexed_AAAAAAAA" = 0L),
        maxUMIs = c("multiplexed_AAAAAAAA" = Inf),
        minGenes = c("multiplexed_AAAAAAAA" = 0L),
        maxGenes = c("multiplexed_AAAAAAAA" = Inf),
        maxMitoRatio = c("multiplexed_AAAAAAAA" = 0L),
        minNovelty = c("multiplexed_AAAAAAAA" = 0L)
    )
    expect_identical(ncol(x), nCells)
})
