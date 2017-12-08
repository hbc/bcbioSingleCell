context("selectSamples")

load(system.file(
    file.path("inst", "extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("inst", "extdata", "filtered.rda"),
    package = "bcbioSingleCell"))

test_that("selectSamples", {
    # Add a quiet argument here in the future
    subset <- suppressMessages(
        selectSamples(filtered, sampleName = "M1_seq_rep_1")
    )
    expect_is(subset, "bcbioSingleCell")
    expect_identical(
        metadata(subset)[["selectSamples"]],
        TRUE
    )
    expect_identical(
        dim(subset),
        c(19525L, 470L)
    )
    expect_identical(
        sampleMetadata(subset)[["sampleID"]],
        factor(c("run1_AGAGGATA"))
    )
    # Ensure that bcbio cellular barcodes get updated correctly
    cb <- bcbio(subset, "cellularBarcodes")
    ids1 <- sampleMetadata(subset) %>%
        .[["sampleID"]] %>%
        as.character()
    ids2 <- names(cb)
    expect_identical(ids1, ids2)
    expect_identical(
        lapply(cb, class) %>%
            unlist() %>%
            unique(),
        "data.frame"
    )
})

test_that("Match failure", {
    expect_equal(
        suppressWarnings(
            selectSamples(filtered, sampleID = "XXX")
        ),
        NULL
    )
    expect_warning(
        selectSamples(filtered, sampleID = "XXX"),
        "Match failure: sampleID = XXX"
    )
    expect_warning(
        selectSamples(filtered, sampleID = "XXX"),
        "No samples matched"
    )
})

test_that("Stop on unfiltered object", {
    expect_error(
        selectSamples(bcb, interestingGroups = "homozygote"),
        "'filterCells\\(\\)' hasn't been applied to this dataset"
    )
})
