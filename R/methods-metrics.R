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
#' @importFrom dplyr left_join mutate_all
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
        if (isTRUE(aggregateReplicates)) {
            .checkAggregate(metadata, stop = TRUE)
            metadata[["sampleName"]] <- metadata[["sampleNameAggregate"]]
            metadata[["sampleNameAggregate"]] <- NULL
        }

        # Define the cell2sample mappings
        # This uses a stashed `data.frame` as of v0.0.22, for better speed
        cell2sample <- metadata(object)[["cell2sample"]]
        if (is.null(cell2sample)) {
            cell2sample <- .cell2sample(
                cells = colData[["cellID"]],
                samples = metadata[["sampleID"]]
            )
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
