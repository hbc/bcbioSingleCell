context("Read Functions")



# readCellRanger ===============================================================
test_that("readCellRanger", {
    object <- suppressWarnings(readCellRanger(
        uploadDir = system.file(
            "extdata/cellranger",
            package = "bcbioSingleCell"
        ),
        organism = "Homo sapiens"
    ))
    expect_is(object, "SingleCellExperiment")
    expect_identical(dim(object), c(500L, 500L))
    expect_identical(
        sampleNames(object),
        c(pbmc_1 = "pbmc")
    )
    expect_identical(
        object = head(colnames(object), n = 1L),
        expected = "pbmc_1_AAACCTGAGAAGGCCT"
    )
    expect_identical(
        object = head(rownames(object), n = 1L),
        expected = "ENSG00000004487"
    )
    expect_identical(
        lapply(metadata(object), class),
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
