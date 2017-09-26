#' Sample Barcode Metrics
#'
#' @rdname metrics
#' @name metrics
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams AllGenerics
#' @param filterCells Only return cells that pass filtering criteria.
#'
#' @return [data.frame] with cellular barcodes as rows.
NULL



# Constructors ====
#' Metrics Constructor
#'
#' @inheritParams AllGenerics
#' @param filterCells Return only metrics for the filtered cells.
#'
#' @return [data.frame].
#' @noRd
.metricsBCB <- function(object) {
    colData <- colData(object) %>%
        as.data.frame() %>%
        rownames_to_column()

    # Get the barcode format based on the umiType
    umiType <- metadata(object)[["umiType"]]
    if (umiType == "surecell") {
        barcode <- "[ACGT]{6}_[ACGT]{6}_[ACGT]{6}"
    } else {
        barcode <- "[ACGT]{8}_[ACGT]{8}"
    }

    # Check for presence of valid barcode in cell identifiers
    if (!all(str_detect(colData[["rowname"]], paste0(barcode, "$")))) {
        stop("Failed to detect surecell barcodes")
    }

    # Match the sampleID from the cellular barcode (rowname)
    match <- colData[["rowname"]] %>%
        str_match(paste0("^(.+)_(", barcode, ")$")) %>%
        as.data.frame() %>%
        set_colnames(c("rowname", "sampleID", "cellularBarcode")) %>%
        # Remove the unnecessary cellularBarcode column
        mutate(cellularBarcode = NULL)

    # In a future update, attempt to detect duplicate columns present in colData
    # and sampleMetadata. We may want to remove the duplicated columns present
    # in the sampleMetadata before joining with the colData.

    # Join the sample metadata by sampleID, extracted from the barcode ID
    left_join(colData, match, by = "rowname") %>%
        left_join(sampleMetadata(object), by = "sampleID") %>%
        as.data.frame() %>%
        column_to_rownames()
}



# Methods ====
#' @rdname metrics
#' @export
setMethod("metrics", "bcbioSingleCell", function(
    object,
    filterCells = TRUE) {
    if (isTRUE(filterCells)) {
        cells <- metadata(object)[["filterCells"]]
        if (!is.null(cells)) {
            object <- object[, cells]
        }
    }
    .metricsBCB(object)
})



#' @rdname metrics
#' @export
setMethod("metrics", "bcbioSingleCellLegacy", .metricsBCB)



#' @rdname metrics
#' @export
setMethod("metrics", "seurat", function(object) {
    object@data.info %>%
        camel()
})
