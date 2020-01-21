#' @name updateObject
#' @author Michael Steinbaugh
#' @note Updated 2020-01-20.
#'
#' @inherit BiocGenerics::updateObject
#' @inheritParams acidroxygen::params
#'
#' @examples
#' data(bcb)
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



## Updated 2020-01-20.
`updateObject,bcbioSingleCell` <-  # nolint
    function(object) {
        cli_h1("Update object")
        sce <- as(object, "SingleCellExperiment")
        cells <- colnames(sce)
        assays <- assays(sce)
        colData <- colData(sce)
        metadata <- metadata(sce)
        version <- metadata[["version"]]
        assert(is(version, c("package_version", "numeric_version")))
        cli_text(sprintf(
            fmt = "Upgrading {.var bcbioSingleCell} from version %s to %s.",
            as.character(version),
            as.character(.version)
        ))

        ## Assays --------------------------------------------------------------
        cli_h2("Assays")
        ## Ensure raw counts are always named "counts".
        if ("assay" %in% names(assays)) {
            ## Versions < 0.1 (e.g. 0.0.21).
            cli_alert("Renaming {.var assay} to {.var counts}.")
            names(assays)[names(assays) == "assay"] <- "counts"
        } else if ("raw" %in% names(assays)) {
            cli_alert("Renaming {.var raw} assay to {.var counts}.")
            names(assays)[names(assays) == "raw"] <- "counts"
        }
        assays <- Filter(Negate(is.null), assays)
        ## Put the required assays first, in order.
        assays <- assays[unique(c(requiredAssays, names(assays)))]
        assert(isSubset(requiredAssays, names(assays)))

        ## Column data metrics -------------------------------------------------
        cli_h2("Column data")
        ## Update legacy column names.
        if (isSubset(c("nCount", "nUMI"), colnames(colData))) {
            cli_alert("Renaming {.var nCount} to {.var nRead}.")
            colnames(colData)[colnames(colData) == "nCount"] <- "nRead"
            cli_alert("Renaming {.var nUMI} to {.var nCount}.")
            colnames(colData)[colnames(colData) == "nUMI"] <- "nCount"
        }
        if (isSubset("nGene", colnames(colData))) {
            cli_alert("Renaming {.var nGene} to {.var nFeature}.")
            colnames(colData)[colnames(colData) == "nGene"] <- "nFeature"
            cli_alert(paste(
                "Renaming {.var log10GenesPerUMI} to",
                "{.var log10FeaturesPerCount}."
            ))
            colnames(colData)[colnames(colData) == "log10GenesPerUMI"] <-
                "log10FeaturesPerCount"
        }
        ## Move sampleData into colData.
        ## Require that all `sampleData` columns are now slotted in `colData`.
        if ("sampleData" %in% names(metadata)) {
            sampleData <- metadata[["sampleData"]]
        } else if ("sampleMetadata" %in% names(metadata)) {
            sampleData <- metadata[["sampleMetadata"]]
        } else {
            sampleData <- NULL
        }
        if (!is.null(sampleData)) {
            cli_alert(paste(
                "Moving {.var sampleData} from {.fun metadata} into",
                "{.fun colData}."
            ))
            assert(isSubset("sampleID", colnames(sampleData)))
            ## Starting using `DataFrame` in place of `data.frame` in v0.1.7.
            sampleData <- as(sampleData, "DataFrame")
            colData <- colData[
                ,
                setdiff(colnames(colData), colnames(sampleData)),
                drop = FALSE
            ]
            cli_alert("Mapping cells to samples.")
            cell2sample <- mapCellsToSamples(
                cells = cells,
                samples = as.character(sampleData[["sampleID"]])
            )
            assert(is.factor(cell2sample))
            colData[["sampleID"]] <- cell2sample
            sampleData[["sampleID"]] <- rownames(sampleData)
            colData <- leftJoin(x = colData, y = sampleData, by = "sampleID")
            assert(
                is(colData, "DataFrame"),
                identical(rownames(colData), colnames(object))
            )
            ## Ensure rows are ordered to match the object.
            colData <- colData[cells, , drop = FALSE]
        }

        ## Metadata ------------------------------------------------------------
        cli_h2("Metadata")
        ## dataVersions
        dataVersions <- metadata[["dataVersions"]]
        if (is(dataVersions, "data.frame")) {
            cli_alert("Setting {.var dataVersions} as {.var DataFrame}.")
            metadata[["dataVersions"]] <- as(dataVersions, "DataFrame")
        }
        ## ensemblRelease
        if ("ensemblVersion" %in% names(metadata)) {
            cli_alert(
                "Renaming {.var ensemblVersion} to {.var ensemblRelease}."
            )
            names(metadata)[
                names(metadata) == "ensemblVersion"] <- "ensemblRelease"
        }
        if (
            is.numeric(metadata[["ensemblRelease"]]) &&
            !is.integer(metadata[["ensemblRelease"]])
        ) {
            cli_alert("Setting {.var ensemblRelease} as integer.")
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
            cli_alert("Renaming {.var gtfFile} to {.var gffFile}.")
            names(metadata)[names(metadata) == "gtfFile"] <- "gffFile"
        }
        if (!"gffFile" %in% names(metadata)) {
            cli_alert("Setting {.var gffFile} as empty character.")
            metadata[["gffFile"]] <- character()
        }
        ## lanes
        if (!is.integer(metadata[["lanes"]])) {
            cli_alert("Setting {.var lanes} as integer.")
            metadata[["lanes"]] <- as.integer(metadata[["lanes"]])
        }
        ## level
        if (!"level" %in% names(metadata)) {
            cli_alert("Setting {.var level} as genes.")
            metadata[["level"]] <- "genes"
        }
        ## programVersions
        if (!"programVersions" %in% names(metadata) &&
            "programs" %in% names(metadata)) {
            cli_alert("Renaming {.var programs} to {.var programVersions}.")
            names(metadata)[names(metadata) == "programs"] <- "programVersions"
        }
        programVersions <- metadata[["programVersions"]]
        if (is(programVersions, "data.frame")) {
            metadata[["programVersions"]] <- as(programVersions, "DataFrame")
        }
        ## sampleMetadataFile
        if (!is.character(metadata[["sampleMetadataFile"]])) {
            cli_alert("Setting {.var sampleMetadataFile} as empty character.")
            metadata[["sampleMetadataFile"]] <- character()
        }
        ## sessionInfo
        ## Support for legacy devtoolsSessionInfo stash.
        ## Previously, we stashed both devtools* and utils* variants.
        if ("devtoolsSessionInfo" %in% names(metadata)) {
            cli_alert("Simplifying stashed {.var sessionInfo}.")
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
        cat_line()
        cli_alert_success(
            "Update of {.var bcbioSingleCell} object was successful."
        )
        bcb
    }



#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature("bcbioSingleCell"),
    definition = `updateObject,bcbioSingleCell`
)
