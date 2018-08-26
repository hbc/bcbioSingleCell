setValidity(
    "bcbioSingleCell",
    function(object) {
        stopifnot(metadata(object)[["version"]] >= 0.1)
        stopifnot(!.hasSlot(object, "bcbio"))

        # Assays ---------------------------------------------------------------
        assert_is_subset("counts", names(assays(object)))

        # Row data -------------------------------------------------------------
        assert_is_all_of(rowRanges(object), "GRanges")
        assert_is_all_of(rowData(object), "DataFrame")

        # Column data ----------------------------------------------------------
        # Require that metrics columns are defined
        assert_is_subset(metricsCols, colnames(colData(object)))

        # Ensure that sample-level metadata is also defined at cell-level.
        # We're doing this in long format in the colData slot.
        assert_is_subset(
            x = colnames(sampleData(object, interestingGroups = NULL)),
            y = colnames(colData(object))
        )

        # Ensure that `interestingGroups` isn't slotted in colData
        stopifnot(!"interestingGroups" %in% colnames(colData(object)))

        # Don't allow the user to manually define `sampleID` column
        stopifnot(!"sampleID" %in% colnames(sampleData(object)))

        # Check that the levels set in `sampleData()` match `colData()`
        sampleDataLevels <- lapply(
            X = sampleData(object, interestingGroups = NULL),
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
            X = colData(object)[, names(sampleDataLevels)],
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
