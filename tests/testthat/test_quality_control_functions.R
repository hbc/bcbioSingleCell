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
        x <- filterCells(
            indrops_small,
            minUMIs = 0L,
            maxUMIs = Inf,
            minGenes = 0L,
            maxGenes = Inf,
            maxMitoRatio = 1L,
            minNovelty = 0L,
            minCellsPerGene = 0L
        )
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
            minUMIs = c(rep_1 = 0L),
            maxUMIs = c(rep_1 = Inf),
            minGenes = c(rep_1 = 0L),
            maxGenes = c(rep_1 = Inf),
            minNovelty = c(rep_1 = 0L),
            maxMitoRatio = c(rep_1 = 0L),
            minCellsPerGene = 10L
        )
    )
})



# metrics ======================================================================
test_that("metrics : bcbioSingleCell", {
    x <- metrics(indrops_small)
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            nCount = "integer",
            nUMI = "integer",
            nGene = "integer",
            nCoding = "integer",
            nMito = "integer",
            log10GenesPerUMI = "numeric",
            mitoRatio = "numeric",
            sampleName = "factor",
            fileName = "factor",
            description = "factor",
            index = "factor",
            sequence = "factor",
            aggregate = "factor",
            revcomp = "factor",
            "interestingGroups" = "factor"
        )
    )
})

test_that("metrics : seurat", {
    x <- metrics(Seurat::pbmc_small)
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            nGene = "integer",
            nUMI = "numeric",
            origIdent = "factor",
            res0x8 = "character",
            res1 = "character",
            sampleName = "factor",
            interestingGroups = "factor",
            ident = "factor"
        )
    )
})
