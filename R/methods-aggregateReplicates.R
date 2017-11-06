#' Aggregate Replicates
#'
#' @rdname aggregateReplicates
#' @name aggregateReplicates
#' @family Data Management Utilities
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom basejump aggregateReplicates
#'
#' @inherit basejump::aggregateReplicates
NULL



# Constructors ====
#' @importFrom dplyr filter mutate select
#' @importFrom parallel mclapply
#' @importFrom stringr str_match
#' @importFrom tibble rownames_to_column
.aggregateReplicates <- function(
    object) {
    cells <- metrics(
        object,
        aggregateReplicates = TRUE,
        filterCells = FALSE) %>%
        select(c("sampleID", "sampleName", "cellularBarcode")) %>%
        mutate(
            cells = paste(
                .data[["sampleName"]],
                .data[["cellularBarcode"]],
                sep = "_"
            )
        ) %>%
        pull("cells")

    # Aggregate the counts
    counts <- aggregateReplicates(assay(object), cells = cells)
    # Check that the count number of counts matches
    if (!identical(sum(assay(object)), sum(counts))) {
        stop("Aggregated counts sum doens't match the original",
             call. = FALSE)
    }

    # Recalculate cellular barcode metrics
    annotable <- annotable(object)
    prefilter <- metadata(object)[["prefilter"]]
    metrics <- calculateMetrics(
        counts,
        annotable = annotable,
        prefilter = prefilter)
    if (isTRUE(prefilter)) {
        # Subset the counts matrix to match the metrics
        counts <- counts[, rownames(metrics)]
    }

    # Update the metadata slot
    metadata <- metadata(object)
    sampleMetadata <- sampleMetadata(object)
    metadata[["sampleMetadata"]] <- sampleMetadata
    metadata[["allSamples"]] <- FALSE

    # Return bcbioSingleCell
    se <- SummarizedExperiment(
        assays = list(assay = counts),
        rowData = annotable,
        colData = metrics,
        metadata = metadata
    )
    new("bcbioSingleCell", se)
}



# Methods ====
#' @rdname aggregateReplicates
#' @export
setMethod(
    "aggregateReplicates",
    signature("bcbioSingleCell"),
    .aggregateReplicates)
