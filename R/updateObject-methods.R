#' @name updateObject
#' @author Michael Steinbaugh
#' @note Updated 2019-07-24.
#'
#' @inherit BiocGenerics::updateObject
#' @inheritParams basejump::params
#'
#' @examples
#' data(indrops)
#' x <- updateObject(indrops)
#' print(x)
NULL



## Updated 2019-07-24.
`updateObject,bcbioSingleCell` <-  # nolint
    function(object, rowRanges = NULL) {
        assert(isAny(rowRanges, classes = c("GRanges", "NULL")))

        assays <- slot(object, name = "assays")
        cells <- colnames(assays[[1L]])
        if (is.null(rowRanges)) {
            rowRanges <- slot(object, name = "rowRanges")
        }
        colData <- slot(object, name = "colData")
        metadata <- slot(object, name = "metadata")

        version <- metadata[["version"]]
        assert(is(version, c("package_version", "numeric_version")))
        message(paste0("Upgrading from ", version, " to ", .version, "."))

        ## Assays --------------------------------------------------------------
        ## Coerce `ShallowSimpleListAssays` S4 class to standard list.
        names <- names(assays)
        assays <- lapply(seq_along(assays), function(a) {
            assay <- assays[[a]]
            assert(identical(colnames(assay), cells))
            assay
        })
        names(assays) <- names
        rm(names)

        ## Ensure raw counts are always named "counts".
        if ("assay" %in% names(assays)) {
            ## Versions < 0.1 (e.g. 0.0.21).
            message("Renaming `assay` to `counts`.")
            assays[["counts"]] <- assays[["assay"]]
            assays[["assay"]] <- NULL
        } else if ("raw" %in% names(assays)) {
            message("Renaming `raw` assay to `counts`.")
            assays[["counts"]] <- assays[["raw"]]
            assays[["raw"]] <- NULL
        }

        assays <- Filter(Negate(is.null), assays)
        ## Put the required assays first, in order
        assays <- assays[unique(c(requiredAssays, names(assays)))]
        assert(isSubset(requiredAssays, names(assays)))

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
            message("Moving `sampleData` columns to `colData`.")
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
                samples = sampleData[["sampleID"]]
            )
            assert(is.factor(cell2sample))
            colData[["rowname"]] <- rownames(colData)
            colData[["sampleID"]] <- cell2sample
            sampleData[["sampleID"]] <- rownames(sampleData)
            colData <- merge(
                x = colData,
                y = sampleData,
                by = "sampleID",
                all.x = TRUE
            )
            rownames(colData) <- colData[["rowname"]]
            colData[["rowname"]] <- NULL
            ## Ensure rows are ordered to match the object.
            colData <- colData[cells, , drop = FALSE]
        }

        ## Metadata ------------------------------------------------------------
        ## dataVersions
        dataVersions <- metadata[["dataVersions"]]
        if (is(dataVersions, "data.frame")) {
            message("Setting dataVersions as DataFrame.")
            metadata[["dataVersions"]] <- as(dataVersions, "DataFrame")
        }

        ## ensemblRelease
        if ("ensemblVersion" %in% names(metadata)) {
            message("Renaming ensemblVersion to ensemblRelease.")
            metadata[["ensemblRelease"]] <- metadata[["ensemblVersion"]]
            metadata[["ensemblVersion"]] <- NULL
        }
        if (
            is.numeric(metadata[["ensemblRelease"]]) &&
            !is.integer(metadata[["ensemblRelease"]])
        ) {
            message("Setting ensemblRelease as integer.")
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
            message("Renaming gtfFile to gffFile.")
            metadata[["gffFile"]] <- metadata[["gtfFile"]]
            metadata[["gtfFile"]] <- NULL
        }
        if (!"gffFile" %in% names(metadata)) {
            message("Setting gffFile as empty character")
            metadata[["gffFile"]] <- character()
        }

        ## lanes
        if (!is.integer(metadata[["lanes"]])) {
            message("Setting lanes as integer.")
            metadata[["lanes"]] <- as.integer(metadata[["lanes"]])
        }

        ## level
        if (!"level" %in% names(metadata)) {
            message("Setting level as genes.")
            metadata[["level"]] <- "genes"
        }

        ## programVersions
        if (!"programVersions" %in% names(metadata) &&
            "programs" %in% names(metadata)) {
            message("Renaming programs to programVersions.")
            metadata[["programVersions"]] <- metadata[["programs"]]
            metadata <- metadata[setdiff(names(metadata), "programs")]
        }
        programVersions <- metadata[["programVersions"]]
        if (is(programVersions, "data.frame")) {
            metadata[["programVersions"]] <- as(programVersions, "DataFrame")
        }

        ## sampleMetadataFile
        if (!is.character(metadata[["sampleMetadataFile"]])) {
            message("Setting sampleMetadataFile as empty character.")
            metadata[["sampleMetadataFile"]] <- character()
        }

        ## Drop legacy slots.
        keep <- setdiff(
            x = names(metadata),
            y = c("cell2sample", "sampleData", "sampleMetadata")
        )
        metadata <- metadata[keep]

        ## Return --------------------------------------------------------------
        .new.bcbioSingleCell(
            assays = assays,
            rowRanges = rowRanges,
            colData = colData,
            metadata = metadata
        )
    }



#' @rdname updateObject
#' @export
setMethod(
    f = "updateObject",
    signature = signature("bcbioSingleCell"),
    definition = `updateObject,bcbioSingleCell`
)
