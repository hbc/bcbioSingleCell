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
#' @importFrom dplyr left_join mutate mutate_all
#' @importFrom magrittr set_colnames
#' @importFrom parallel mclapply
#' @importFrom stringr str_match
#' @importFrom tibble column_to_rownames rownames_to_column
#' @export
setMethod(
    "metrics",
    signature("bcbioSingleCell"),
    function(
        object,
        interestingGroups,
        filterCells = FALSE,
        aggregateReplicates = FALSE) {
        if (isTRUE(filterCells)) {
            object <- .applyFilterCutoffs(object)
        }

        colData <- colData(object) %>%
            as.data.frame() %>%
            rownames_to_column("cellID")

        # Prepare the metadata
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
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

        # Check for existence of `sampleID` in the `cellID`
        sampleID <- metadata[["sampleID"]]
        cellID <- colData[["cellID"]]
        cell2sample <- mclapply(seq_along(sampleID), function(a) {
            pattern <- paste0("^(", sampleID[[a]], ")_([ACGT_]{6,})(_[0-9]+)?$")
            x <- str_match(cellID, pattern = pattern) %>%
                as.data.frame() %>%
                .[, 1:2] %>%
                set_colnames(c("cellID", "sampleID")) %>%
                filter(!is.na(.data[["sampleID"]]))
        }) %>%
            bind_rows() %>%
            mutate(sampleID = as.factor(.data[["sampleID"]]))
        if (!identical(length(cellID), nrow(cell2sample))) {
            stop("Failed to correctly match sample IDs to cellular barcodes",
                 call. = FALSE)
        }

        colData %>%
            left_join(cell2sample, by = "cellID") %>%
            left_join(metadata, by = "sampleID") %>%
            column_to_rownames("cellID")
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
        interestingGroups <- basejump::interestingGroups(object)
    }
    metrics <- slot(object, "meta.data") %>%
        as.data.frame() %>%
        camel(strict = FALSE) %>%
        rownames_to_column("cellID") %>%
        mutate_if(is.character, as.factor)
    metadata <- sampleMetadata(object)
    suppressMessages(left_join(metrics, metadata)) %>%
        uniteInterestingGroups(interestingGroups) %>%
        column_to_rownames("cellID")
})
