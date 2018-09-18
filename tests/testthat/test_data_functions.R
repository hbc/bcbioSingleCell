context("Data Functions")



# aggregateCols ==========================================================
test_that("aggregateCols", {
    x <- aggregateCols(indrops_small)
    expect_s4_class(x, "SingleCellExperiment")
    expect_identical(
        dim(x),
        c(500L, 500L)
    )
    map <- metadata(x)[["aggregateCols"]]
    expect_is(map, "factor")
    expect_identical(length(map), 500L)
    expect_identical(length(levels(map)), 500L)
})



# gene2symbol ==================================================================
test_that("gene2symbol : bcbioSingleCell", {
    x <- gene2symbol(indrops_small)
    expect_is(x, "gene2symbol")
})



# interestingGroups ============================================================
test_that("interestingGroups", {
    expect_identical(
        object = interestingGroups(indrops_small),
        expected = "sampleName"
    )
})

test_that("interestingGroups<-", {
    expect_silent(
        interestingGroups(indrops_small) <- "sampleName"
    )
    expect_error(
        interestingGroups(indrops_small) <- "XXX"
    )
})



# sampleData ===================================================================
test_that("sampleData", {
    object <- sampleData(indrops_small)
    expect_identical(
        object = lapply(object, class),
        expected = list(
            sampleName  = "factor",
            fileName  = "factor",
            description  = "factor",
            index = "factor",
            sequence = "factor",
            aggregate = "factor",
            revcomp = "factor",
            interestingGroups = "factor"
        )
    )
})



# selectSamples ================================================================
test_that("selectSamples", {
    x <- selectSamples(indrops_small, sampleName = "rep_1")
    expect_s4_class(x, "bcbioSingleCell")
    expect_true(metadata(x)[["selectSamples"]])
    expect_identical(dim(x), c(500L, 500L))
    expect_identical(
        rownames(sampleData(x)),
        "multiplexed_AAAAAAAA"
    )
})

test_that("selectSamples : Match failure", {
    expect_error(
        selectSamples(indrops_small, sampleName = "XXX"),
        "\"sampleName\" metadata column doesn't contain XXX"
    )
})



# subsetPerSample ==============================================================
test_that("subsetPerSample", {
    x <- subsetPerSample(indrops_small, assignAndSave = FALSE)
    expect_is(x, "list")
    expect_identical(names(x), "multiplexed_AAAAAAAA")

    # Assign and save.
    unlink("subsetPerSample", recursive = TRUE)
    subsetPerSample(
        object = indrops_small,
        assignAndSave = TRUE,
        dir = "subsetPerSample"
    )
    expect_identical(
        list.files("subsetPerSample"),
        "multiplexed_AAAAAAAA.rda"
    )
    load("subsetPerSample/multiplexed_AAAAAAAA.rda")
    expect_identical(
        dim(multiplexed_AAAAAAAA),
        c(500L, 500L)
    )
    unlink("subsetPerSample", recursive = TRUE)
})



# topBarcodes ==================================================================
test_that("topBarcodes", {
    # tibble
    x <- topBarcodes(indrops_small, return = "tibble")
    expect_identical(dplyr::group_vars(x), "sampleID")
    expect_identical(
        lapply(x, class),
        list(
            sampleID = "factor",
            sampleName = "factor",
            nUMI = "integer",
            cellID = "character"
        )
    )

    # list
    x <- topBarcodes(indrops_small, return = "list")
    expect_is(x, "list")
})
