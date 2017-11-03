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
#' @param interestingGroups Interesting group, to use for colors.
#' @param filterCells Show only the cells that have passed filtering cutoffs.
#' @param aggregateReplicates Aggregate technical replicates, if specified. This
#'   function uses values assigned in the `sampleNameAggregate` column of the
#'   internal sample metadata [data.frame].
#'
#' @seealso [sampleMetadata()].
#'
#' @return [data.frame] with cellular barcodes as rows.
NULL



# Methods ====
#' @rdname metrics
#' @importFrom basejump uniteInterestingGroups
#' @importFrom dplyr left_join mutate
#' @importFrom magrittr set_colnames
#' @importFrom stringr str_match
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
setMethod(
    "metrics",
    signature("bcbioSingleCell"),
    function(
        object,
        interestingGroups,
        filterCells = TRUE,
        aggregateReplicates = TRUE) {
        # Filter cells that passed QC checks, if desired
        if (isTRUE(filterCells)) {
            cells <- metadata(object)[["filterCells"]]
            if (!is.null(cells)) {
                object <- object[, cells]
            }
        }

        colData <- colData(object) %>%
            as.data.frame() %>%
            rownames_to_column("cellID")

        # Prepare the metadata
        if (missing(interestingGroups)) {
            interestingGroups <- metadata(object)[["interestingGroups"]]
        }
        metadata <- metadata(object)[["sampleMetadata"]] %>%
            as.data.frame() %>%
            uniteInterestingGroups(interestingGroups) %>%
            # Ensure all the columns are factors here. This will also sanitize
            # any previously saved datasets, where the sample identifiers are
            # character vectors.
            mutate_all(as.factor)

        # Rename `sampleName` when aggregating replicates
        if (isTRUE(aggregateReplicates) &
            "sampleNameAggregate" %in% colnames(metadata)) {
            metadata[["sampleName"]] <- metadata[["sampleNameAggregate"]]
            metadata[["sampleNameAggregate"]] <- NULL
        }

        # Join the sample metadata by maching against the cellular barcode
        # Get the barcode format based on `umiType`
        umiType <- metadata(object)[["umiType"]]
        if (umiType == "surecell") {
            barcode <- "[ACGT]{6}_[ACGT]{6}_[ACGT]{6}"
        } else {
            barcode <- "[ACGT]{8}_[ACGT]{8}"
        }
        # Check for presence of valid barcode in cell identifiers
        if (!all(grepl(x = colData[["cellID"]],
                       pattern = paste0(barcode, "$")))) {
            stop("Failed to detect surecell barcodes", call. = FALSE)
        }
        # Match `sampleID` from the cellular barcode (rowname)
        match <- colData[["cellID"]] %>%
            str_match(paste0("^(.+)_(", barcode, ")$")) %>%
            as.data.frame() %>%
            set_colnames(c("cellID", "sampleID", "cellularBarcode")) %>%
            # Remove the unnecessary cellularBarcode column
            mutate(
                cellularBarcode = NULL,
                sampleID = as.factor(.data[["sampleID"]])
            )
        # Join the sample metadata by `sampleID`, extracted from `cellID`
        metrics <- colData %>%
            left_join(match, by = "cellID") %>%
            left_join(metadata, by = "sampleID") %>%
            as.data.frame() %>%
            column_to_rownames("cellID")
        metrics
    })



#' @rdname metrics
#' @importFrom basejump camel
#' @importFrom dplyr mutate_if
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
setMethod(
    "metrics",
    signature("seurat"),
    function(
        object,
        interestingGroups) {
    if (missing(interestingGroups)) {
        # Check to see if bcbio interestingGroups are stashed.
        # Stop if they're not defined.
        interestingGroups <- slot(object, "misc") %>%
            .[["bcbio"]] %>%
            .[["interestingGroups"]]
        if (is.null(interestingGroups)) {
            stop("'interestingGroups' is required", call. = FALSE)
        }
    }
    # This slot was previously named `data.info` in earlier releases
    metrics <- slot(object, "meta.data") %>%
        as.data.frame() %>%
        rownames_to_column() %>%
        camel(strict = FALSE)
    # Check for required columns, and make sure they're factors
    requiredCols <- c("description", "sampleID", "sampleName")
    if (!all(requiredCols %in% colnames(metrics))) {
        stop(paste(
            "object@meta.data must contain these columns:",
            toString(requiredCols)
        ), call. = FALSE)
    }
    # Ensure all strings are factors
    metrics <- mutate_if(metrics, is.character, as.factor)
    # Create the `interestingGroups` column required for QC plots
    metrics <- uniteInterestingGroups(metrics, interestingGroups)
    metrics <- column_to_rownames(metrics)
    metrics
})
