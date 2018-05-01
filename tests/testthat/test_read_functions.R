context("Read Functions")



# readCellRanger ===============================================================
test_that("readCellRanger", {
    extdataDir <- system.file("extdata", package = "bcbioSingleCell")
    uploadDir <- file.path(extdataDir, "cellranger")
    refdataDir <- file.path(uploadDir, "refdata-cellranger-hg19-1.2.0")
    sampleMetadataFile <- file.path(uploadDir, "metadata.csv")
    x <- readCellRanger(
        uploadDir = uploadDir,
        refdataDir = refdataDir,
        sampleMetadataFile = sampleMetadataFile
    )
    expect_is(x, "SingleCellExperiment")
    expect_identical(dim(x), c(500L, 500L))
    expect_identical(
        colnames(x)[1L:2L],
        c("aggregation_AAACCTGGTTTACTCT_1",
          "aggregation_AAACGGGGTATCTGCA_1")
    )
    expect_identical(
        rownames(x)[1L:2L],
        c("ENSG00000005022",
          "ENSG00000008018")
    )
    expect_identical(
        lapply(metadata(x), class),
        list(
            "version" = c("package_version", "numeric_version"),
            "pipeline" = "character",
            "level" = "character",
            "uploadDir" = "character",
            "sampleDirs" = "character",
            "sampleMetadataFile" = "character",
            "sampleData" = "data.frame",
            "interestingGroups" = "character",
            "cell2sample" = "factor",
            "organism" = "character",
            "genomeBuild" = "character",
            "ensemblRelease" = "integer",
            "rowRangesMetadata" = c("tbl_df", "tbl", "data.frame"),
            "umiType" = "character",
            "allSamples" = "logical",
            "refdataDir" = "character",
            "refJSON" = "list",
            "call" = "call",
            "date" = "Date",
            "wd" = "character",
            "utilsSessionInfo" = "sessionInfo",
            "devtoolsSessionInfo" = "session_info"
        )
    )
})



# readCellTypeMarkers ==========================================================
test_that("readCellTypeMarkers : Mus musculus", {
    file <- system.file(
        "extdata/cell_type_markers.csv",
        package = "bcbioSingleCell"
    )
    gene2symbol <- makeGene2symbolFromEnsembl("Homo sapiens")
    x <- readCellTypeMarkers(
        file = file,
        gene2symbol = gene2symbol
    )
    group <- dplyr::group_vars(x)
    expect_identical(group, "cellType")
    x <- x[1L, ]
    y <- tibble::tibble(
        "cellType" = "B Cell",
        "geneID" = "ENSG00000177455",
        "geneName" = "CD19"
    ) %>%
        dplyr::group_by(!!rlang::sym("cellType"))
    expect_identical(x, y)
})
