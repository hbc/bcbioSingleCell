#' bcbio single-cell RNA-seq data set
#'
#' `bcbioSingleCell` is an S4 class that extends `SingleCellExperiment`, and is
#' designed to store a bcbio single-cell RNA-seq analysis. This class contains
#' read counts saved as a sparse matrix (`sparseMatrix`), sample metadata, and
#' cell quality control metrics.
#'
#' @author Michael Steinbaugh, Rory Kirchner
#' @note Updated 2019-10-30.
#' @export
setClass(
    Class = "bcbioSingleCell",
    contains = "SingleCellExperiment"
)
setValidity(
    Class = "bcbioSingleCell",
    method = function(object) {
        colData <- colData(object)
        metadata <- metadata(object)
        sampleData <- sampleData(object)
        ## Return invalid for all objects older than v0.1.
        version <- metadata[["version"]]
        ok <- validate(
            is(version, "package_version"),
            version >= 0.1
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ## Check for legacy bcbio slot.
        ok <- validate(!.hasSlot(object, "bcbio"))
        if (!isTRUE(ok)) {
            return(ok)
        }
        ## Assays --------------------------------------------------------------
        ok <- validate(isSubset("counts", names(assays(object))))
        if (!isTRUE(ok)) {
            return(ok)
        }
        ## Row data ------------------------------------------------------------
        ok <- validate(
            is(rowRanges(object), "GenomicRanges"),
            is(rowData(object), "DataFrame")
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ## Column data ---------------------------------------------------------
        ok <- validate(
            ## Require that metrics columns are defined.
            isSubset(.metricsCols, colnames(colData)),
            ## Ensure that `interestingGroups` isn't slotted in colData.
            areDisjointSets("interestingGroups", colnames(colData))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ## Metadata ------------------------------------------------------------
        df <- c("DFrame", "DataFrame")
        ok <- validateClasses(
            object = metadata,
            expected = list(
                allSamples = "logical",
                bcbioCommandsLog = "character",
                bcbioLog = "character",
                dataVersions = df,
                date = "Date",
                ensemblRelease = "integer",
                genomeBuild = "character",
                gffFile = "character",
                interestingGroups = "character",
                lanes = "integer",
                level = "character",
                organism = "character",
                pipeline = "character",
                programVersions = df,
                projectDir = "character",
                runDate = "Date",
                sampleDirs = "character",
                sampleMetadataFile = "character",
                sessionInfo = "session_info",
                umiType = "character",
                uploadDir = "character",
                version = "package_version",
                wd = "character",
                yaml = "list"
            ),
            subset = TRUE
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        ## Check that level is defined.
        ok <- validate(
            !isSubset("sampleName", names(metadata)),
            isSubset(metadata[["level"]], c("genes", "transcripts"))
        )
        if (!isTRUE(ok)) {
            return(ok)
        }
        TRUE
    }
)
