context("Read Functions")



# readCellRanger ===============================================================
test_that("readCellRanger", {
    x <- suppressWarnings(readCellRanger(
        uploadDir = system.file(
            "extdata/cellranger",
            package = "bcbioSingleCell"
        ),
        organism = "Homo sapiens"
    ))
    expect_is(x, "SingleCellExperiment")
    expect_identical(dim(x), c(500L, 500L))
    expect_identical(
        sampleNames(x),
        c(pbmc_1 = "pbmc")
    )
    expect_identical(
        colnames(x)[[1L]],
        "pbmc_1_AAACCTGAGAAGGCCT"
    )
    expect_identical(
        rownames(x)[[1L]],
        "ENSG00000243485"
    )
    expect_identical(
        lapply(metadata(x), class),
        list(
            version = c("package_version", "numeric_version"),
            pipeline = "character",
            level = "character",
            uploadDir = "character",
            sampleDirs = "character",
            sampleMetadataFile = "character",
            sampleData = "data.frame",
            interestingGroups = "character",
            cell2sample = "factor",
            organism = "character",
            genomeBuild = "character",
            ensemblRelease = "integer",
            rowRangesMetadata = c("tbl_df", "tbl", "data.frame"),
            umiType = "character",
            allSamples = "logical",
            call = "call",
            date = "Date",
            wd = "character",
            utilsSessionInfo = "sessionInfo",
            devtoolsSessionInfo = "session_info"
        )
    )
})
