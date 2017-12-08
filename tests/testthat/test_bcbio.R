context("bcbio")

load(system.file(
    file.path("extdata", "bcb.rda"),
    package = "bcbioSingleCell"))
load(system.file(
    file.path("extdata", "seurat.rda"),
    package = "bcbioSingleCell"))

test_that("bcbioSingleCell", {
    slot <- bcbio(bcb)
    expect_is(slot, "SimpleList")
    expect_identical(
        lapply(slot, class),
        list(
            "cellularBarcodes" = "list"
        )
    )
})

test_that("bcbioSingleCell assignment", {
    data <- bcb
    bcbio(data, "stash") <- "test"
    expect_identical(
        bcbio(data, "stash"),
        "test"
    )
})

test_that("seurat", {
    tibble <- c("tbl_df", "tbl", "data.frame")
    slot <- bcbio(seurat)
    expect_is(slot, "list")
    expect_identical(
        lapply(slot, class),
        list("version" = c("package_version", "numeric_version"),
             "pipeline" = "character",
             "uploadDir" = "character",
             "sampleDirs" = "character",
             "sampleMetadataFile" = "character",
             "sampleMetadata" = "data.frame",
             "interestingGroups" = "character",
             "cell2sample" = "factor",
             "organism" = "character",
             "genomeBuild" = "character",
             "ensemblVersion" = "NULL",
             "gtfFile" = "NULL",
             "annotable" = "data.frame",
             "gene2symbol" = "data.frame",
             "umiType" = "character",
             "allSamples" = "logical",
             "prefilter" = "logical",
             "projectDir" = "character",
             "template" = "character",
             "runDate" = "Date",
             "yamlFile" = "character",
             "yaml" = "list",
             "tx2gene" = "logical",
             "dataVersions" = tibble,
             "programs" = tibble,
             "bcbioLog" = "character",
             "bcbioCommandsLog" = "character",
             "cellularBarcodeCutoff" = "numeric",
             "date" = "Date",
             "wd" = "character",
             "utilsSessionInfo" = "sessionInfo",
             "devtoolsSessionInfo" = "session_info",
             "aggregateReplicates" = "factor",
             "filterCells" = "character",
             "filterGenes" = "character",
             "filterParams" = "numeric",
             "subset" = "logical")
    )
})

test_that("seurat assignment", {
    data <- seurat
    bcbio(data, "stash") <- "test"
    expect_identical(
        bcbio(data, "stash"),
        "test"
    )
})
