#' bcbio single-cell RNA-seq data set
#'
#' `bcbioSingleCell` is an S4 class that extends `SingleCellExperiment`, and is
#' designed to store a bcbio single-cell RNA-seq analysis. This class contains
#' read counts saved as a sparse matrix (`sparseMatrix`), sample metadata, and
#' cell quality control metrics.
#'
#' @family S4 classes
#' @author Michael Steinbaugh
#' @export
#'
#' @seealso `bcbioSingleCell`.
setClass(
    Class = "bcbioSingleCell",
    contains = "SingleCellExperiment"
)
setValidity(
    Class = "bcbioSingleCell",
    method = function(object) {
        valid <- list()

        colData <- colData(object)
        metadata <- metadata(object)
        sampleData <- sampleData(object)

        # Return invalid for all objects older than v0.1.
        version <- metadata[["version"]]
        valid[["version"]] <- validate(
            is(version, "package_version"),
            version >= 0.1
        )

        valid[["general"]] <- validate(
            !.hasSlot(object, "bcbio")
        )

        # Assays ---------------------------------------------------------------
        valid[["assays"]] <- validate(
            isSubset("counts", names(assays(object)))
        )

        # Row data -------------------------------------------------------------
        valid[["rowData"]] <- validate(
            is(rowRanges(object), "GRanges"),
            is(rowData(object), "DataFrame")
        )

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

        valid[["colData"]] <- validate(
            # Require that metrics columns are defined.
            isSubset(metricsCols, colnames(colData)),
            # Ensure that `interestingGroups` isn't slotted in colData.
            areDisjointSets("interestingGroups", colnames(colData)),
            # Ensure that sample-level metadata is also defined at cell-level.
            # We're doing this in long format in the colData slot.
            isSubset(colnames(sampleData), colnames(colData)),
            identical(sampleDataLevels, colDataLevels)
        )

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
        valid[["metadata"]] <- validateClasses(
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

        valid[["metadata2"]] <- validate(
            !isSubset("sampleName", names(metadata)),
            isSubset(metadata[["level"]], c("genes", "transcripts"))
        )

        .valid(list = valid)
    }
)



.valid <- function(list) {
    invalid <- Filter(f = Negate(isTRUE), x = list)
    if (hasLength(invalid)) {
        unlist(invalid)
    } else {
        TRUE
    }
}
