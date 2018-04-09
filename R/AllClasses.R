setOldClass(Classes = c("grouped_df", "tbl_df", "tibble"))



#' bcbioSingleCell Class
#'
#' `bcbioSingleCell` extends `SingleCellExperiment` and designed to store a
#' bcbio single-cell RNA-seq analysis. This class contains read counts saved as
#' a sparse matrix (`dgCMatrix`), sample metadata, and cell quality control
#' metrics.
#'
#' @note `bcbioSingleCell` extended `SummarizedExperiment` prior to v0.1.0,
#'   where we migrated to `SingleCellExperiment`.
#'
#' @author Michael Steinbaugh
#'
#' @seealso
#' - [loadSingleCell()], [loadCellRanger()].
#' - [SingleCellExperiment::SingleCellExperiment()].
#' - `.S4methods(class = "bcbioSingleCell")`.
#'
#' @export
bcbioSingleCell <- setClass(
    "bcbioSingleCell",
    contains = "SingleCellExperiment"
)



# Validity =====================================================================
setValidity(
    "bcbioSingleCell",
    function(object) {
        stopifnot(metadata(object)[["version"]] >= 0.1)
        assert_is_all_of(object, "SingleCellExperiment")
        stopifnot(!.hasSlot(object, "bcbio"))

        # Assays ===============================================================
        assert_are_identical("counts", names(assays(object)))

        # Row data =============================================================
        assert_is_all_of(rowRanges(object), "GRanges")
        assert_is_all_of(rowData(object), "DataFrame")
        # Require gene-to-symbol mappings
        assert_is_subset(
            x = c("geneID", "geneName"),
            y = colnames(rowData(object))
        )

        # Column data ==========================================================
        # Check that all of the columns are numeric
        colDataCheck <- vapply(
            X = slot(object, "colData"),
            FUN = is.numeric,
            FUN.VALUE = logical(1L),
            USE.NAMES = TRUE
        )
        if (!all(colDataCheck)) {
            stop(paste(
                paste(
                    "Non-numeric colData columns:",
                    toString(names(colDataCheck[!colDataCheck]))
                ),
                bcbioBase::updateMessage,
                sep = "\n"
            ))
        }

        # Metadata =============================================================
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
        # - dataVersions: tbl_df
        # - gffFile: character
        # - programVersions: tbl_df
        # - projectDir: character
        # - runDate: Date
        # - template: character
        # - yaml: list
        #
        # v0.2.0
        # - loadCellRanger: call
        # - loadSingleCell: call
        # - txdb: TxDb

        # Class checks
        requiredMetadata <- list(
            "allSamples" = "logical",
            "cell2sample" = "factor",
            "date" = "Date",
            "devtoolsSessionInfo" = "session_info",
            "ensemblRelease" = "integer",
            "genomeBuild" = "character",
            "interestingGroups" = "character",
            "isSpike" = "character",
            "level" = "character",
            "organism" = "character",
            "pipeline" = "character",
            "sampleData" = "data.frame",
            "sampleDirs" = "character",
            "sampleMetadataFile" = "character",
            "umiType" = "character",
            "uploadDir" = "character",
            "utilsSessionInfo" = "sessionInfo",
            "version" = "package_version",
            "wd" = "character"
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
            print(classChecks)
            stop(paste(
                "Metadata class checks failed.",
                bcbioBase::updateMessage,
                sep = "\n"
            ))
        }

        # level
        assert_is_subset(
            x = metadata[["level"]],
            y = c("genes", "transcripts")
        )

        # sampleData
        invisible(lapply(metadata[["sampleData"]], assert_is_factor))

        TRUE
    }
)
