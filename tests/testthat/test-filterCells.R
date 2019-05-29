context("filterCells")

test_that("No filtering", {
    # Expecting an object with the same dimensions by default.
    invisible(capture.output(
        x <- filterCells(indrops)
    ))
    expect_s4_class(x, "bcbioSingleCell")
    expect_identical(dim(x), dim(indrops))
})

test_that("Expected cutoff failure", {
    expect_error(
        filterCells(indrops, minUMIs = Inf),
        "No cells passed `minUMIs` cutoff"
    )
})

with_parameters_test_that(
    "Parameterized cutoff tests", {
        args[["object"]] <- indrops
        invisible(capture.output(
            x <- do.call(what = filterCells, args = args)
        ))
        expect_s4_class(x, "bcbioSingleCell")
        expect_is(metadata(x)[["filterParams"]], "list")
        expect_is(metadata(x)[["filterCells"]], "character")
        expect_is(metadata(x)[["filterGenes"]], "character")
        expect_identical(metadata(x)[["subset"]], TRUE)
        expect_identical(dim(x), dim)
    },
    # Refer to the quality control R Markdown for actual recommended cutoffs.
    # These are skewed, and designed to work with our minimal dataset.
    args = list(
        list(minUMIs = 2000L),
        list(maxUMIs = 2500L),
        list(minGenes = 45L),
        list(maxGenes = 49L),
        list(maxMitoRatio = 0.1),
        list(minNovelty = 0.5),
        list(minCellsPerGene = 95L)
    ),
    dim = list(
        c(50L, 35L),
        c(50L, 88L),
        c(50L, 95L),
        c(50L, 81L),
        c(50L, 22L),
        c(50L, 81L),
        c(45L, 100L)
    )
)

test_that("Per sample cutoffs", {
    # Get the count of sample1 (run1_AGAGGATA)
    # We're applying no filtering to that sample
    sampleNames <- sampleNames(indrops)
    expect_identical(
        sampleNames,
        c(multiplexed_AAAAAAAA = "rep_1")
    )
    invisible(capture.output(
        object <- filterCells(
            object = indrops,
            minUMIs = c(rep_1 = 1L),
            maxUMIs = c(rep_1 = Inf),
            minGenes = c(rep_1 = 1L),
            maxGenes = c(rep_1 = Inf),
            maxMitoRatio = c(rep_1 = 1L),
            minNovelty = c(rep_1 = 0L)
        )
    ))
    expect_identical(
        object = metadata(object)[["filterParams"]],
        expected = list(
            nCells = Inf,
            minUMIs = c(rep_1 = 1L),
            maxUMIs = c(rep_1 = Inf),
            minGenes = c(rep_1 = 1L),
            maxGenes = c(rep_1 = Inf),
            minNovelty = c(rep_1 = 0L),
            maxMitoRatio = c(rep_1 = 1L),
            minCellsPerGene = 1L
        )
    )
})
