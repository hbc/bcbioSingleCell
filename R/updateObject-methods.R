#' Update object
#'
#' @name updateObject
#' @author Michael Steinbaugh
#' @note Updated 2022-05-09.
#'
#' @inheritParams AcidRoxygen::params
#'
#' @return Modified object.
#'
#' @examples
#' data(bcb)
#'
#' ## bcbioSingleCell ====
#' updateObject(bcb)
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



## Updated 2022-05-09.
`updateObject,bcbioSingleCell` <- # nolint
    function(object, ..., verbose = FALSE) {
        assert(isFlag(verbose))
        if (isTRUE(verbose)) {
            h1("Update object")
        }
        sce <- as(object, "SingleCellExperiment")
        cells <- colnames(sce)
        assays <- assays(sce)
        rowRanges <- rowRanges(sce)
        colData <- colData(sce)
        metadata <- metadata(sce)
        version <- metadata[["version"]]
        assert(is(version, "package_version"))
        if (isTRUE(verbose)) {
            alert(sprintf(
                fmt = "Upgrading {.var %s} from version %s to %s.",
                "bcbioSingleCell",
                as.character(version),
                as.character(.version)
            ))
        }
        ## Assays --------------------------------------------------------------
        if (isTRUE(verbose)) {
            h2("Assays")
        }
        ## Ensure raw counts are always named "counts".
        if (isSubset("assay", names(assays))) {
            ## Versions < 0.1 (e.g. 0.0.21).
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "assay", "counts"
                ))
            }
            names(assays)[names(assays) == "assay"] <- "counts"
        } else if (isSubset("raw", names(assays))) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} assay to {.var %s}.",
                    "raw", "counts"
                ))
            }
            names(assays)[names(assays) == "raw"] <- "counts"
        }
        assays <- Filter(Negate(is.null), assays)
        ## Put the required assays first, in order.
        assays <- assays[unique(c(.requiredAssays, names(assays)))]
        assert(isSubset(.requiredAssays, names(assays)))
        ## Row data ------------------------------------------------------------
        if (hasNames(mcols(rowRanges))) {
            mcols(rowRanges) <-
                camelCase(mcols(rowRanges), strict = TRUE)
        }
        ## Column data ---------------------------------------------------------
        if (isTRUE(verbose)) {
            h2("Column data")
        }
        colnames(colData) <- camelCase(colnames(colData), strict = TRUE)
        if (isSubset(c("nCount", "nUmi"), colnames(colData))) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "nCount", "nRead"
                ))
            }
            colnames(colData)[colnames(colData) == "nCount"] <- "nRead"
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "nUmi", "nCount"
                ))
            }
            colnames(colData)[colnames(colData) == "nUmi"] <- "nCount"
        }
        if (isSubset("nGene", colnames(colData))) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "nGene", "nFeature"
                ))
            }
            colnames(colData)[colnames(colData) == "nGene"] <- "nFeature"
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "log10GenesPerUmi", "log10FeaturesPerCount"
                ))
            }
            colnames(colData)[colnames(colData) == "log10GenesPerUmi"] <-
                "log10FeaturesPerCount"
        }
        ## Move sampleData into colData.
        if (isSubset("sampleData", names(metadata))) {
            sampleData <- metadata[["sampleData"]]
        } else if (isSubset("sampleMetadata", names(metadata))) {
            sampleData <- metadata[["sampleMetadata"]]
        } else {
            sampleData <- NULL
        }
        if (!is.null(sampleData)) {
            colnames(sampleData) <-
                camelCase(colnames(sampleData), strict = TRUE)
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Moving {.var %s} from {.fun %s} into {.fun %s}.",
                    "sampleData", "metadata", "colData"
                ))
            }
            assert(isSubset("sampleId", colnames(sampleData)))
            sampleData <- as(sampleData, "DataFrame")
            colData <- colData[
                ,
                setdiff(colnames(colData), colnames(sampleData)),
                drop = FALSE
            ]
            if (isTRUE(verbose)) {
                alert("Mapping cells to samples.")
            }
            cell2sample <- mapCellsToSamples(
                cells = cells,
                samples = as.character(sampleData[["sampleId"]])
            )
            assert(is.factor(cell2sample))
            colData[["sampleId"]] <- cell2sample
            sampleData[["sampleId"]] <- rownames(sampleData)
            colData <- leftJoin(x = colData, y = sampleData, by = "sampleId")
            assert(
                is(colData, "DataFrame"),
                identical(rownames(colData), colnames(object))
            )
            ## Ensure rows are ordered to match the object.
            colData <- colData[cells, , drop = FALSE]
        }
        ## Metadata ------------------------------------------------------------
        if (isTRUE(verbose)) {
            h2("Metadata")
        }
        ## dataVersions.
        dataVersions <- metadata[["dataVersions"]]
        if (is(dataVersions, "data.frame")) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} as {.cls %s}.",
                    "dataVersions", "DataFrame"
                ))
            }
            metadata[["dataVersions"]] <- as(dataVersions, "DataFrame")
        }
        ## ensemblRelease.
        if (isSubset("ensemblVersion", names(metadata))) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "ensemblVersion", "ensemblRelease"
                ))
            }
            names(metadata)[
                names(metadata) == "ensemblVersion"
            ] <- "ensemblRelease"
        }
        if (
            is.numeric(metadata[["ensemblRelease"]]) &&
                !is.integer(metadata[["ensemblRelease"]])
        ) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} as integer.",
                    "ensemblRelease"
                ))
            }
            metadata[["ensemblRelease"]] <-
                as.integer(metadata[["ensemblRelease"]])
        }
        ## Update the version, if necessary.
        if (!identical(metadata[["version"]], .version)) {
            metadata[["originalVersion"]] <- metadata[["version"]]
            metadata[["version"]] <- .version
        }
        ## gffFile.
        if (isSubset("gtfFile", names(metadata))) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "gtfFile", "gffFile"
                ))
            }
            names(metadata)[names(metadata) == "gtfFile"] <- "gffFile"
        }
        if (!isSubset("gffFile", names(metadata))) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "gffFile", "empty character"
                ))
            }
            metadata[["gffFile"]] <- character()
        }
        ## lanes.
        if (!is.integer(metadata[["lanes"]])) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "lanes", "integer"
                ))
            }
            metadata[["lanes"]] <- as.integer(metadata[["lanes"]])
        }
        ## level.
        if (!isSubset("level", names(metadata))) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "level", "genes"
                ))
            }
            metadata[["level"]] <- "genes"
        }
        ## programVersions.
        if (!isSubset("programVersions", names(metadata)) &&
            isSubset("programs", names(metadata))
        ) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Renaming {.var %s} to {.var %s}.",
                    "programs", "programVersions"
                ))
            }
            names(metadata)[names(metadata) == "programs"] <- "programVersions"
        }
        programVersions <- metadata[["programVersions"]]
        if (is(programVersions, "data.frame")) {
            metadata[["programVersions"]] <- as(programVersions, "DataFrame")
        }
        ## sampleMetadataFile.
        if (!is.character(metadata[["sampleMetadataFile"]])) {
            if (isTRUE(verbose)) {
                alert(sprintf(
                    "Setting {.var %s} as {.val %s}.",
                    "sampleMetadataFile", "empty character"
                ))
            }
            metadata[["sampleMetadataFile"]] <- character()
        }
        ## sessionInfo.
        if (isSubset("utilsSessionInfo", names(metadata))) {
            if (isTRUE(verbose)) {
                alert(sprintf("Simplifying stashed {.var %s}.", "sessionInfo"))
            }
            names(metadata)[
                names(metadata) == "utilsSessionInfo"
            ] <- "sessionInfo"
            metadata[["devtoolsSessionInfo"]] <- NULL
        }
        ## Drop legacy slots.
        keep <- setdiff(
            x = names(metadata),
            y = c("cell2sample", "sampleData", "sampleMetadata")
        )
        metadata <- metadata[keep]
        ## Return --------------------------------------------------------------
        assays(sce) <- assays
        rowRanges(sce) <- rowRanges
        colData(sce) <- colData
        metadata(sce) <- metadata
        bcb <- new(Class = "bcbioSingleCell", sce)
        validObject(bcb)
        if (isTRUE(verbose)) {
            alertSuccess(sprintf(
                "Update of {.var %s} object was successful.",
                "bcbioSingleCell"
            ))
        }
        bcb
    }



#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature(object = "bcbioSingleCell"),
    definition = `updateObject,bcbioSingleCell`
)
