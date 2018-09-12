# tibble
setOldClass(Classes = c("grouped_df", "tbl_df"))



# bcbioSingleCell ==============================================================
#' `bcbioSingleCell` Class
#'
#' `bcbioSingleCell` is an S4 class that extends `SingleCellExperiment`, and is
#' designed to store a bcbio single-cell RNA-seq analysis. This class contains
#' read counts saved as a sparse matrix (`dgCMatrix`), sample metadata, and cell
#' quality control metrics.
#'
#' @export
setClass(
    Class = "bcbioSingleCell",
    contains = "SingleCellExperiment"
)

setValidity(
    Class = "bcbioSingleCell",
    method = function(object) {
        stopifnot(metadata(object)[["version"]] >= 0.1)
        stopifnot(!.hasSlot(object, "bcbio"))

        # Assays ---------------------------------------------------------------
        assert_is_subset("counts", names(assays(object)))

        # Row data -------------------------------------------------------------
        assert_is_all_of(rowRanges(object), "GRanges")
        assert_is_all_of(rowData(object), "DataFrame")

        # Column data ----------------------------------------------------------
        colData <- colData(object)
        sampleData <- sampleData(object)
        sampleData[["interestingGroups"]] <- NULL

        # Require that metrics columns are defined.
        assert_is_subset(metricsCols, colnames(colData))

        # Ensure that `interestingGroups` isn't slotted in colData.
        assert_are_disjoint_sets("interestingGroups", colnames(colData))

        # Ensure that sample-level metadata is also defined at cell-level.
        # We're doing this in long format in the colData slot.
        assert_is_subset(colnames(sampleData), colnames(colData))

        # Check that the levels set in `sampleData()` match `colData()`
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
            X = colData[, names(sampleDataLevels)],
            FUN = levels
        )
        assert_are_identical(sampleDataLevels, colDataLevels)

        # Metadata -------------------------------------------------------------
        metadata <- metadata(object)

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

        # Legacy metadata:
        # Switched to DataFrame in v0.1.17.
        # Now using `colData()` directly in v0.2.2.
        # sampleData = c("DataFrame", "data.frame")

        # Class checks
        requiredMetadata <- list(
            allSamples = "logical",
            cell2sample = "factor",
            date = "Date",
            devtoolsSessionInfo = "session_info",
            ensemblRelease = "integer",
            genomeBuild = "character",
            interestingGroups = "character",
            level = "character",
            organism = "character",
            pipeline = "character",
            sampleDirs = "character",
            sampleMetadataFile = "character",
            umiType = "character",
            uploadDir = "character",
            utilsSessionInfo = "sessionInfo",
            version = "package_version",
            wd = "character"
        )
        classChecks <- invisible(mapply(
            name <- names(requiredMetadata),
            expected <- requiredMetadata,
            MoreArgs = list(metadata = metadata),
            FUN = function(name, expected, metadata) {
                actual <- class(metadata[[name]])
                if (!length(intersect(expected, actual))) {
                    FALSE
                } else {
                    TRUE
                }
            },
            SIMPLIFY = TRUE,
            USE.NAMES = TRUE
        ))
        if (!all(classChecks)) {
            stop(paste(
                "Metadata class checks failed.",
                updateMessage,
                printString(classChecks),
                sep = "\n"
            ))
        }

        # level
        assert_is_subset(
            x = metadata[["level"]],
            y = c("genes", "transcripts")
        )

        TRUE
    }
)



# CellRanger ===================================================================
#' `CellRanger` Class
#'
#' Extends `SingleCellExperiment`, with additional validity checks on the
#' `metadata()` slot.
#'
#' @export
setClass(
    Class = "CellRanger",
    contains = "SingleCellExperiment"
)

setValidity(
    Class = "CellRanger",
    method = function(object) {
        TRUE
    }
)
