#' @name updateObject
#' @author Michael Steinbaugh
#' @note Updated 2019-08-12.
#'
#' @inherit BiocGenerics::updateObject
#' @inheritParams acidroxygen::params
#'
#' @examples
#' data(indrops)
#' updateObject(indrops)
#'
#' ## Example that depends on remote file.
#' ## > x <- import(
#' ## >     file = file.path(
#' ## >         bcbioSingleCellTestsURL,
#' ## >         "bcbioSingleCell_0.1.0.rds"
#' ## >     )
#' ## > )
#' ## > x <- updateObject(x)
#' ## > x
NULL



## Updated 2019-08-12.
`updateObject,bcbioSingleCell` <-  # nolint
    function(object) {
        sce <- as(object, "SingleCellExperiment")

        cells <- colnames(sce)
        assays <- assays(sce)
        colData <- colData(sce)
        metadata <- metadata(sce)

        version <- metadata[["version"]]
        assert(is(version, c("package_version", "numeric_version")))
        message(sprintf(
            fmt = "Upgrading bcbioSingleCell from version %s to %s.",
            as.character(version),
            as.character(.version)
        ))

        ## Assays --------------------------------------------------------------
        ## Ensure raw counts are always named "counts".
        if ("assay" %in% names(assays)) {
            ## Versions < 0.1 (e.g. 0.0.21).
            message("Renaming 'assay' to 'counts'.")
            names(assays)[names(assays) == "assay"] <- "counts"
        } else if ("raw" %in% names(assays)) {
            message("Renaming 'raw' assay to 'counts'.")
            names(assays)[names(assays) == "raw"] <- "counts"
        }

        assays <- Filter(Negate(is.null), assays)
        ## Put the required assays first, in order.
        assays <- assays[unique(c(requiredAssays, names(assays)))]
        assert(isSubset(requiredAssays, names(assays)))

        ## Column data metrics -------------------------------------------------
        ## Update legacy column names.
        if (isSubset(c("nCount", "nUMI"), colnames(colData))) {
            message("Renaming 'nCount' to 'nRead'.")
            colnames(colData)[colnames(colData) == "nCount"] <- "nRead"
            message("Renaming 'nUMI' to 'nCount'.")
            colnames(colData)[colnames(colData) == "nUMI"] <- "nCount"
        }
        if (isSubset("nGene", colnames(colData))) {
            message("Renaming 'nGene' to 'nFeature'.")
            colnames(colData)[colnames(colData) == "nGene"] <- "nFeature"
            message("Renaming 'log10GenesPerUMI' to 'log10FeaturesPerCount'.")
            colnames(colData)[colnames(colData) == "log10GenesPerUMI"] <-
                "log10FeaturesPerCount"
        }

        ## Move sampleData into colData ----------------------------------------
        ## Require that all `sampleData` columns are now slotted in `colData`.
        if ("sampleData" %in% names(metadata)) {
            sampleData <- metadata[["sampleData"]]
        } else if ("sampleMetadata" %in% names(metadata)) {
            sampleData <- metadata[["sampleMetadata"]]
        } else {
            sampleData <- NULL
        }
        if (!is.null(sampleData)) {
            message("Moving 'sampleData' from 'metadata()' into 'colData()'.")
            assert(isSubset("sampleID", colnames(sampleData)))
            ## Starting using `DataFrame` in place of `data.frame` in v0.1.7.
            sampleData <- as(sampleData, "DataFrame")
            colData <- colData[
                ,
                setdiff(colnames(colData), colnames(sampleData)),
                drop = FALSE
            ]
            message("Mapping cells to samples.")
            cell2sample <- mapCellsToSamples(
                cells = cells,
                samples = as.character(sampleData[["sampleID"]])
            )
            assert(is.factor(cell2sample))
            colData[["sampleID"]] <- cell2sample
            sampleData[["sampleID"]] <- rownames(sampleData)
            colData <- left_join(x = colData, y = sampleData, by = "sampleID")
            assert(
                is(colData, "DataFrame"),
                identical(rownames(colData), colnames(object))
            )
            ## Ensure rows are ordered to match the object.
            colData <- colData[cells, , drop = FALSE]
        }

        ## Metadata ------------------------------------------------------------
        ## dataVersions
        dataVersions <- metadata[["dataVersions"]]
        if (is(dataVersions, "data.frame")) {
            message("Setting 'dataVersions' as DataFrame.")
            metadata[["dataVersions"]] <- as(dataVersions, "DataFrame")
        }

        ## ensemblRelease
        if ("ensemblVersion" %in% names(metadata)) {
            message("Renaming 'ensemblVersion' to 'ensemblRelease'.")
            names(metadata)[
                names(metadata) == "ensemblVersion"] <- "ensemblRelease"
        }
        if (
            is.numeric(metadata[["ensemblRelease"]]) &&
            !is.integer(metadata[["ensemblRelease"]])
        ) {
            message("Setting 'ensemblRelease' as integer.")
            metadata[["ensemblRelease"]] <-
                as.integer(metadata[["ensemblRelease"]])
        }

        ## Update the version, if necessary.
        if (!identical(metadata[["version"]], .version)) {
            metadata[["originalVersion"]] <- metadata[["version"]]
            metadata[["version"]] <- .version
        }

        ## gffFile
        if ("gtfFile" %in% names(metadata)) {
            message("Renaming 'gtfFile' to 'gffFile'.")
            names(metadata)[names(metadata) == "gtfFile"] <- "gffFile"
        }
        if (!"gffFile" %in% names(metadata)) {
            message("Setting 'gffFile' as empty character.")
            metadata[["gffFile"]] <- character()
        }

        ## lanes
        if (!is.integer(metadata[["lanes"]])) {
            message("Setting 'lanes' as integer.")
            metadata[["lanes"]] <- as.integer(metadata[["lanes"]])
        }

        ## level
        if (!"level" %in% names(metadata)) {
            message("Setting 'level' as genes.")
            metadata[["level"]] <- "genes"
        }

        ## programVersions
        if (!"programVersions" %in% names(metadata) &&
            "programs" %in% names(metadata)) {
            message("Renaming 'programs' to 'programVersions'.")
            names(metadata)[names(metadata) == "programs"] <- "programVersions"
        }
        programVersions <- metadata[["programVersions"]]
        if (is(programVersions, "data.frame")) {
            metadata[["programVersions"]] <- as(programVersions, "DataFrame")
        }

        ## sampleMetadataFile
        if (!is.character(metadata[["sampleMetadataFile"]])) {
            message("Setting 'sampleMetadataFile' as empty character.")
            metadata[["sampleMetadataFile"]] <- character()
        }

        ## sessionInfo
        ## Support for legacy devtoolsSessionInfo stash.
        ## Previously, we stashed both devtools* and utils* variants.
        if ("devtoolsSessionInfo" %in% names(metadata)) {
            message("Simplifying stashed 'sessionInfo'.")
            names(metadata)[
                names(metadata) == "devtoolsSessionInfo"] <- "sessionInfo"
            metadata[["utilsSessionInfo"]] <- NULL
        }

        ## Drop legacy slots.
        keep <- setdiff(
            x = names(metadata),
            y = c("cell2sample", "sampleData", "sampleMetadata")
        )
        metadata <- metadata[keep]

        ## Return --------------------------------------------------------------
        ## Note that this approach will keep `spikeNames` set, if defined.
        assays(sce) <- assays
        colData(sce) <- colData
        metadata(sce) <- metadata
        bcb <- new(Class = "bcbioSingleCell", sce)
        validObject(bcb)
        bcb
    }



#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature("bcbioSingleCell"),
    definition = `updateObject,bcbioSingleCell`
)
