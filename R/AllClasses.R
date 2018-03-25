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
#' @slot bcbio `SimpleList` containing additional bcbio run data with dimensions
#' that don't match the count matrix. This is currently used to store all
#' unfiltered cellular barcodes for quality control analysis.
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
        assert_has_dimnames(object)
        stopifnot(!.hasSlot(object, "bcbio"))

        # Assays ===============================================================
        assert_are_identical("raw", names(assays(object)))

        # Row data =============================================================
        assert_is_all_of(rowRanges(object), "GRanges")
        assert_is_all_of(rowData(object), "data.frame")
        # Require gene-to-symbol mappings
        assert_is_subset(
            x = c("geneID", "geneName"),
            y = colnames(rowData(object))
        )

        # Column data ==========================================================
        # Check that all of the columns are numeric
        colDataCheck <- vapply(
            X = colData(object),
            FUN = is.numeric,
            FUN.VALUE = logical(1L),
            USE.NAMES = TRUE
        )
        if (!all(colDataCheck)) {
            abort(paste(
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
        # - lanes
        # - tx2gene
        #
        # v0.2.0
        # - loadCellRanger: call
        # - loadSingleCell: call
        # - txdb: TxDb

        # Class checks (order independent)
        requiredMetadata <- list(
            "allSamples" = "logical",
            "bcbioCommandsLog" = "character",
            "bcbioLog" = "character",
            "cell2sample" = "factor",
            "dataVersions" = "tbl_df",
            "date" = "Date",
            "devtoolsSessionInfo" = "session_info",
            "ensemblRelease" = "integer",
            "genomeBuild" = "character",
            "gffFile" = "character",
            "interestingGroups" = "character",
            "isSpike" = "character",
            "level" = "character",
            "organism" = "character",
            "pipeline" = "character",
            "programVersions" = "tbl_df",
            "projectDir" = "character",
            "rowRangesMetadata" = "tbl_df",
            "runDate" = "Date",
            "sampleData" = "data.frame",
            "sampleDirs" = "character",
            "sampleMetadataFile" = "character",
            "template" = "character",
            "umiType" = "character",
            "unannotatedRows" = "character",
            "uploadDir" = "character",
            "utilsSessionInfo" = "sessionInfo",
            "version" = "package_version",
            "wd" = "character",
            "yaml" = "list"
        )
        classChecks <- invisible(vapply(
            X = seq_along(requiredMetadata),
            FUN = function(a) {
                name <- names(requiredMetadata)[[a]]
                actual <- class(metadata[[name]])
                expected <- requiredMetadata[[a]]
                if (!length(intersect(expected, actual))) {
                    warn(paste(
                        name, "is not", toString(expected)
                    ))
                    FALSE
                } else {
                    TRUE
                }
            },
            FUN.VALUE = logical(1L),
            USE.NAMES = FALSE
        ))
        if (!all(classChecks)) {
            abort(paste(
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
# object@bcbio$cellularBarcodes needs to go in metadata
