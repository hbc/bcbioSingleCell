#' Update an Object to Its Current Class Definition
#'
#' @name updateObject
#' @family S4 Object
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics updateObject
#' @export
#'
#' @inherit BiocGenerics::updateObject
#'
#' @return Valid object.
#'
#' @examples
#' x <- updateObject(indrops_small)
#' print(x)
NULL



#' @rdname updateObject
#' @export
setMethod(
    "updateObject",
    signature("bcbioSingleCell"),
    function(object) {
        version <- slot(object, "metadata")[["version"]]
        assert_is_all_of(version, c("package_version", "numeric_version"))
        message(paste("Upgrading from", version, "to", packageVersion))

        # Coerce to SingleCellExperiment
        sce <- as(object, "SingleCellExperiment")
        validObject(sce)

        assays <- slot(sce, "assays")
        rowRanges <- rowRanges(sce)
        colData <- colData(sce)
        metadata <- metadata(sce)

        # Assays ---------------------------------------------------------------
        # Coerce ShallowSimpleListAssays to list
        assays <- lapply(seq_along(assays), function(a) {
            assay <- assays[[a]]
            assert_are_identical(colnames(assay), colnames(sce))
            assay
        })
        names(assays) <- assayNames(sce)

        # raw counts
        if ("raw" %in% names(assays)) {
            message("Renaming `raw` assay to `counts`")
            assays[["counts"]] <- assays[["raw"]]
            assays[["raw"]] <- NULL
        }

        assays <- Filter(Negate(is.null), assays)
        # Put the required assays first, in order
        assays <- assays[unique(c(requiredAssays, names(assays)))]
        assert_is_subset(requiredAssays, names(assays))

        # Column data ----------------------------------------------------------
        # Require that all `sampleData` columns are now slotted in `colData`.
        sampleData <- metadata[["sampleData"]]
        if (!is.null(sampleData)) {
            message("Moving `sampleData()` columns to `colData()`")
            # Starting using `DataFrame` in place of `data.frame` in v0.1.7.
            sampleData <- as(sampleData, "DataFrame")
            colData <- colData[
                ,
                setdiff(colnames(colData), colnames(sampleData)),
                drop = FALSE
                ]
            cell2sample <- cell2sample(sce)
            assert_is_factor(cell2sample)
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
            # Ensure rows are ordered to match the object.
            colData <- colData[colnames(sce), , drop = FALSE]
        }

        # Metadata -------------------------------------------------------------
        metadata[["sampleData"]] <- NULL
        metadata[["previousVersion"]] <- metadata[["version"]]
        metadata[["version"]] <- packageVersion

        # Return ---------------------------------------------------------------
        .new.bcbioSingleCell(
            assays = assays,
            rowRanges = rowRanges(sce),
            colData = colData,
            metadata = metadata
        )
    }
)
