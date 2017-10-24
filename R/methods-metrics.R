#' Sample Barcode Metrics
#'
#' @rdname metrics
#' @name metrics
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom basejump metrics
#'
#' @inheritParams AllGenerics
#'
#' @param filterCells Show only the cells that have passed filtering cutoffs.
#' @param aggregateReplicates Aggregate technical replicates, if specified. This
#'   function uses values assigned in the `sampleNameAggregate` column of the
#'   internal sample metadata [data.frame].
#'
#' @seealso [sampleMetadata()].
#'
#' @return [data.frame] with cellular barcodes as rows.
NULL



# Constructors ====
#' Metrics Constructor
#'
#' @importFrom dplyr left_join mutate
#' @importFrom magrittr set_colnames
#' @importFrom stringr str_match
#' @importFrom tibble column_to_rownames rownames_to_column
#'
#' @inheritParams AllGenerics
#'
#' @param filterCells Return only metrics for the filtered cells.
#'
#' @return [data.frame].
#' @noRd
.metrics <- function(object, aggregateReplicates = TRUE) {
    colData <- colData(object) %>%
        as.data.frame() %>%
        rownames_to_column("cellID")
    meta <- metadata(object)[["sampleMetadata"]] %>%
        as.data.frame()

    # Rename `sampleName` when aggregating replicates
    if (isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(meta)) {
        meta[["sampleName"]] <- meta[["sampleNameAggregate"]]
        meta[["sampleNameAggregate"]] <- NULL
    }

    # Get the barcode format based on the umiType
    umiType <- metadata(object)[["umiType"]]
    if (umiType == "surecell") {
        barcode <- "[ACGT]{6}_[ACGT]{6}_[ACGT]{6}"
    } else {
        barcode <- "[ACGT]{8}_[ACGT]{8}"
    }

    # Check for presence of valid barcode in cell identifiers
    if (!all(grepl(x = colData[["cellID"]], pattern = paste0(barcode, "$")))) {
        stop("Failed to detect surecell barcodes", call. = FALSE)
    }

    # Match the sampleID from the cellular barcode (rowname)
    match <- colData[["cellID"]] %>%
        str_match(paste0("^(.+)_(", barcode, ")$")) %>%
        as.data.frame() %>%
        set_colnames(c("cellID", "sampleID", "cellularBarcode")) %>%
        # Remove the unnecessary cellularBarcode column
        mutate(cellularBarcode = NULL)

    # Join the sample metadata by sampleID, extracted from the barcode ID
    colData %>%
        left_join(match, by = "cellID") %>%
        left_join(meta, by = "sampleID") %>%
        as.data.frame() %>%
        column_to_rownames("cellID")
}



# Methods ====
#' @rdname metrics
#' @export
setMethod("metrics", "bcbioSingleCell", function(
    object,
    filterCells = TRUE,
    aggregateReplicates = TRUE) {
    if (isTRUE(filterCells)) {
        cells <- metadata(object)[["filterCells"]]
        if (!is.null(cells)) {
            object <- object[, cells]
        }
    }
    .metrics(object, aggregateReplicates = aggregateReplicates)
})



#' @rdname metrics
#' @importFrom basejump camel
#' @export
setMethod("metrics", "seurat", function(object) {
    object@data.info %>%
        camel(strict = FALSE) %>%
        .metrics(aggregateReplicates = FALSE)
})
