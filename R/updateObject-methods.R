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



# Methods ======================================================================
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

        # Assays ===============================================================
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

        # Column data ==========================================================
        # Ensure that all `sampleData` columns are now slotted in `colData`
        sampleData <- metadata[["sampleData"]]
        stopifnot(!is.null(sampleData))
        # Starting using `DataFrame` in place of `data.frame` in v0.1.7
        sampleData <- as(sampleData, "DataFrame")
        colData <- colData[
            ,
            setdiff(colnames(colData), colnames(sampleData)),
            drop = FALSE
        ]
        cell2sample <- cell2sample(sce)
        assert_is_factor(cell2sample)
        colData[["cellID"]] <- rownames(colData)
        colData[["sampleID"]] <- cell2sample
        sampleData[["sampleID"]] <- rownames(sampleData)
        colData <- merge(
            x = colData,
            y = sampleData,
            by = "sampleID",
            all.x = TRUE
        )
        rownames(colData) <- colData[["cellID"]]
        # Ensure rows are reordered to match the original object
        colData <- colData[colnames(sce), ]
        colData[["cellID"]] <- NULL
        sampleData[["sampleID"]] <- NULL

        # Metadata =============================================================
        metadata[["sampleData"]] <- sampleData
        metadata[["previousVersion"]] <- metadata[["version"]]
        metadata[["version"]] <- packageVersion

        # Return ===============================================================
        .new.bcbioSingleCell(
            assays = assays,
            rowRanges = rowRanges(sce),
            colData = colData,
            metadata = metadata
        )
    }
)
