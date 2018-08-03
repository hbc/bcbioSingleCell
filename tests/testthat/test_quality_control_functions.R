context("Quality Control Functions")



# filterCells ==================================================================
test_that("filterCells", {
    invisible(capture.output(
        x <- filterCells(
            object = indrops_small,
            minUMIs = 1000L,
            maxUMIs = Inf,
            minGenes = 100L,
            maxGenes = Inf,
            maxMitoRatio = 0.1,
            minNovelty = 0.7,
            minCellsPerGene = 3L
        )
    ))
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(dim(x), c(500L, 205L))
    expect_is(metadata(x)[["filterParams"]], "list")
    expect_is(metadata(x)[["filterCells"]], "character")
    expect_is(metadata(x)[["filterGenes"]], "character")
    expect_identical(metadata(x)[["subset"]], TRUE)
})

test_that("filterCells : Maximum parameters", {
    # This should return an object with the same dimensions
    invisible(capture.output(
        x <- filterCells(indrops_small)
    ))
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(dim(x), dim(indrops_small))
})

test_that("filterCells : Cutoff failures", {
    expect_error(
        filterCells(indrops_small, minUMIs = Inf),
        "No cells passed `minUMIs` cutoff"
    )
})

test_that("filterCells : Per sample cutoffs", {
    # Get the count of sample1 (run1_AGAGGATA)
    # We're applying no filtering to that sample
    sampleNames <- sampleNames(indrops_small)
    expect_identical(
        sampleNames,
        c(multiplexed_AAAAAAAA = "rep_1")
    )
    invisible(capture.output(
        x <- filterCells(
            object = indrops_small,
            minUMIs = c(rep_1 = 0L),
            maxUMIs = c(rep_1 = Inf),
            minGenes = c(rep_1 = 0L),
            maxGenes = c(rep_1 = Inf),
            maxMitoRatio = c(rep_1 = 0L),
            minNovelty = c(rep_1 = 0L)
        )
    ))
    expect_identical(
        metadata(x)[["filterParams"]],
        list(
            nCells = Inf,
            minUMIs = c(rep_1 = 0L),
            maxUMIs = c(rep_1 = Inf),
            minGenes = c(rep_1 = 0L),
            maxGenes = c(rep_1 = Inf),
            minNovelty = c(rep_1 = 0L),
            maxMitoRatio = c(rep_1 = 0L),
            minCellsPerGene = 1L
        )
    )
})

test_that("filterCells: zinbwave mode", {
    x <- filterCells(indrops_small, zinbwave = TRUE)
    expect_s4_class(x, "bcbioSingleCell")
    expect_true(
        all(c("counts", "normalizedValues", "weights") %in% assayNames(x))
    )
})



# metrics ======================================================================
test_that("metrics", {
    object <- indrops_small
    x <- metrics(object)
    expect_true(all(
        colnames(colData(object)) %in% colnames(x)
    ))
    expect_true("interestingGroups" %in% colnames(x))
})
