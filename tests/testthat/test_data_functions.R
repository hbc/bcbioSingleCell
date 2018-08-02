context("Data Functions")



# aggregateReplicates ==========================================================
test_that("aggregateReplicates", {
    x <- aggregateReplicates(indrops_small)
    expect_s4_class(x, "SingleCellExperiment")
    expect_identical(
        dim(x),
        c(500L, 500L)
    )
    map <- metadata(x)[["aggregateReplicates"]]
    expect_is(map, "factor")
    expect_identical(length(map), 500L)
    expect_identical(length(levels(map)), 500L)
})



# gene2symbol ==================================================================
colnames <- c("geneID", "geneName")

test_that("gene2symbol : bcbioSingleCell", {
    x <- gene2symbol(indrops_small)
    expect_is(x, "data.frame")
    expect_identical(colnames(x), colnames)
})



# interestingGroups ============================================================
test_that("interestingGroups", {
    expect_identical(
        interestingGroups(indrops_small),
        "sampleName"
    )
})

test_that("interestingGroups<-", {
    expect_silent(
        interestingGroups(indrops_small) <- "sampleName"
    )
    expect_error(
        interestingGroups(indrops_small) <- "XXX",
        "The interesting groups \"XXX\" are not defined"
    )
})



# sampleData ===================================================================
all <- list(
    sampleName  = "factor",
    fileName  = "factor",
    description  = "factor",
    index = "factor",
    sequence = "factor",
    aggregate = "factor",
    revcomp = "factor",
    interestingGroups = "factor"
)

test_that("sampleData", {
    x <- sampleData(indrops_small)
    expect_identical(lapply(x, class), all)
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
    # data.frame
    x <- topBarcodes(indrops_small, return = "data.frame")
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
