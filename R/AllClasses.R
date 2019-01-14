#' bcbio single-cell RNA-seq data set
#'
#' `bcbioSingleCell` is an S4 class that extends `SingleCellExperiment`, and is
#' designed to store a bcbio single-cell RNA-seq analysis. This class contains
#' read counts saved as a sparse matrix (`sparseMatrix`), sample metadata, and
#' cell quality control metrics.
#'
#' @author Michael Steinbaugh, Rory Kirchner
#' @export
setClass(
    Class = "bcbioSingleCell",
    contains = "SingleCellExperiment",
    validity = function(object) {
        colData <- colData(object)
        metadata <- metadata(object)
        sampleData <- sampleData(object)

        # Return invalid for all objects older than v0.1.
        version <- metadata[["version"]]
        ok <- validate(
            is(version, "package_version"),
            version >= 0.1
        )
        if (!isTRUE(ok)) return(ok)

        # Check for legacy bcbio slot.
        ok <- validate(!.hasSlot(object, "bcbio"))
        if (!isTRUE(ok)) return(ok)

        # Assays ---------------------------------------------------------------
        ok <- validate(isSubset("counts", names(assays(object))))
        if (!isTRUE(ok)) return(ok)

        # Row data -------------------------------------------------------------
        ok <- validate(
            is(rowRanges(object), "GRanges"),
            is(rowData(object), "DataFrame")
        )
        if (!isTRUE(ok)) return(ok)

        # Column data ----------------------------------------------------------
        sampleData[["interestingGroups"]] <- NULL

        # Check that the levels set in `sampleData` match `colData`
        sampleDataLevels <- lapply(
            X = sampleData,
            FUN = function(x) {
                if (is.factor(x)) {
                    levels(x)
                } else {
                    NULL
                }
            }
        )
        sampleDataLevels <- Filter(Negate(is.null), sampleDataLevels)
        colDataLevels <- lapply(
            X = colData[, names(sampleDataLevels), drop = FALSE],
            FUN = levels
        )

        ok <- validate(
            # Require that metrics columns are defined.
            isSubset(metricsCols, colnames(colData)),
            # Ensure that `interestingGroups` isn't slotted in colData.
            areDisjointSets("interestingGroups", colnames(colData)),
            # Ensure that sample-level metadata is also defined at cell-level.
            # We're doing this in long format in the colData slot.
            isSubset(colnames(sampleData), colnames(colData)),
            identical(sampleDataLevels, colDataLevels)
        )
        if (!isTRUE(ok)) return(ok)

        # Metadata -------------------------------------------------------------
        # Optional metadata:
        # - filterCells
        # - filterGenes
        # - filterParams
        # - filterSummary
        # - lanes: integer
        # - rowRangesMetadata: tbl_df
        # - tx2gene: data.frame
        #
        # bcbio-specific:
        # - bcbioCommandsLog: character
        # - bcbioLog: character
        # - cellularBarcodes: list
        # - dataVersions: tbl_df
        # - gffFile: character
        # - programVersions: tbl_df
        # - projectDir: character
        # - runDate: Date
        # - template: character
        # - yaml: list
        ok <- validateClasses(
            object = metadata,
            expected = list(
                allSamples = "logical",
                date = "Date",
                ensemblRelease = "integer",
                genomeBuild = "character",
                interestingGroups = "character",
                level = "character",
                organism = "character",
                pipeline = "character",
                sampleDirs = "character",
                sampleMetadataFile = "character",
                sessionInfo = "session_info",
                umiType = "character",
                uploadDir = "character",
                version = "package_version",
                wd = "character"
            ),
            subset = TRUE
        )
        if (!isTRUE(ok)) return(ok)

        ok <- validate(
            !isSubset("sampleName", names(metadata)),
            isSubset(metadata[["level"]], c("genes", "transcripts"))
        )
        if (!isTRUE(ok)) return(ok)

        TRUE
    }
)
