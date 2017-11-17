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
#' @inheritParams AllGenerics
NULL



# Constructors ====
#' @importFrom dplyr filter mutate select
#' @importFrom parallel mclapply
#' @importFrom stringr str_match
#' @importFrom tibble rownames_to_column
.aggregateReplicates <- function(object) {
    metadata <- sampleMetadata(object)
    if (!"sampleNameAggregate" %in% colnames(metadata)) {
        stop("'sampleNameAggregate' not present in sample metadata")
    }
    metadata <- metadata[, c("sampleID", "sampleName", "sampleNameAggregate")]
    metadata[["sampleIDAggregate"]] <- make.names(
        metadata[["sampleNameAggregate"]],
        unique = FALSE)
    # We'll end up replacing `sampleID` and `sampleName` columns with the
    # corresponding `*Aggregate` columns.

    metrics <- metrics(object)
    cell2sample <- cell2sample(object)

    map <- left_join(cell2sample, metadata, by = "sampleID")
    cells <- mapply(
        FUN = gsub,
        x = map[["cellID"]],
        pattern = paste0("^", map[["sampleID"]]),
        replacement = map[["sampleIDAggregate"]]
    )

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
    sampleMetadata <- sampleMetadata(object, aggregateReplicates = TRUE)
    metadata[["sampleMetadata"]] <- sampleMetadata

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
