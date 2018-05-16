#' Update Object
#'
#' Update old objects created by the bcbioSingleCell package. The session
#' information metadata is preserved from the time when the bcbio data was
#' originally loaded into R.
#'
#' @name updateObject
#' @family S4 Class Definition
#' @author Michael Steinbaugh
#'
#' @importFrom BiocGenerics updateObject
#'
#' @inheritParams general
#'
#' @return `bcbioRNASeq`.
#'
#' @examples
#' loadRemoteData("http://bcbiosinglecell.seq.cloud/v0.1.6/indrops_small.rda")
#' updateObject(bcb_small)
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

        # Assays ===============================================================
        assays <- slot(sce, "assays")
        if (!all(assayNames(sce) %in% requiredAssays)) {
            # Coerce ShallowSimpleListAssays to list
            assays <- lapply(seq_along(assays), function(a) {
                assay <- assays[[a]]
                if (!identical(colnames(assay), colnames(sce))) {
                    message("Fixing assay colnames")
                    colnames(assay) <- colnames(sce)
                    assay
                }
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
        }

        # Column data ==========================================================
        # Ensure that all `sampleData` columns are now slotted in `colData`
        colData <- colData(sce)
        sampleData <- metadata(sce)[["sampleData"]]
        assert_is_data.frame(sampleData)
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
        colData[["cellID"]] <- NULL
        sampleData[["sampleID"]] <- NULL

        # Metadata =============================================================
        metadata <- metadata(sce)

        # isSpike
        # TODO Add support for updating GRanges using emptyRanges

        # unannotatedRows
        # TODO Check current conventions and update here if necessary

        # version
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
