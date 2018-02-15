context("selectSamples")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "filtered.rda"),
    package = "bcbioSingleCell"))

test_that("selectSamples", {
    # Add a quiet argument here in the future
    subset <- suppressMessages(
        selectSamples(filtered, sampleName = "M1")
    )
    expect_is(subset, "bcbioSingleCell")
    expect_identical(
        metadata(subset)[["selectSamples"]],
        TRUE
    )
    expect_identical(
        dim(subset),
        c(1000L, 243L)
    )
    expect_identical(
        sampleMetadata(subset)[["sampleID"]],
        factor("M1")
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

# FIXME Use assertive for error
test_that("Match failure", {
    expect_error(
        selectSamples(filtered, sampleID = "XXX"),
        "sampleID metadata column doesn't contain XXX"
    )
})

test_that("Stop on unfiltered object", {
    expect_error(
        selectSamples(bcb, interestingGroups = "homozygote"),
        "`filterCells\\(\\)` hasn't been applied to this dataset"
    )
})
