context("loadCellRanger")

test_that("Homo sapiens", {
    extdataDir <- system.file("extdata", package = "bcbioSingleCell")
    uploadDir <- file.path(extdataDir, "cellranger")
    refdataDir <- file.path(extdataDir, "refdata-cellranger-hg19-1.2.0")
    sampleMetadataFile <- file.path(extdataDir, "cellranger.csv")
    # This will warn about missing reference GTF
    bcb <- loadCellRanger(
        uploadDir = uploadDir,
        refdataDir = refdataDir,
        sampleMetadataFile = sampleMetadataFile)
    expect_identical(
        dim(bcb),
        c(500L, 500L)
    )
    expect_identical(
        colnames(bcb)[1L:2L],
        c("aggregation_AAACCTGGTTTACTCT_1",
          "aggregation_AAACGGGGTATCTGCA_1")
    )
    expect_identical(
        rownames(bcb)[1L:2L],
        c("ENSG00000005022",
          "ENSG00000008018")
    )
    expect_identical(
        lapply(metadata(bcb), class),
        list(
            "version" = c("package_version", "numeric_version"),
            "pipeline" = "character",
            "uploadDir" = "character",
            "sampleDirs" = "character",
            "sampleMetadataFile" = "character",
            "sampleMetadata" = "data.frame",
            "interestingGroups" = "character",
            "cell2sample" = "factor",
            "organism" = "character",
            "genomeBuild" = "character",
            "ensemblVersion" = "integer",
            "annotable" = "data.frame",
            "gtfFile" = "character",
            "gene2symbol" = "data.frame",
            "umiType" = "character",
            "allSamples" = "logical",
            "prefilter" = "logical",
            "refdataDir" = "character",
            "refJSON" = "list",
            "date" = "Date",
            "wd" = "character",
            "utilsSessionInfo" = "sessionInfo",
            "devtoolsSessionInfo" = "session_info"
        )
    )
})
