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
#' @importFrom dplyr arrange group_by filter mutate mutate_all select ungroup
#' @importFrom parallel mclapply
#' @importFrom rlang !!! syms
#' @importFrom stats reorder
#' @importFrom stringr str_match
#' @importFrom tibble rownames_to_column
.aggregateReplicates <- function(object) {
    sampleMetadata <- sampleMetadata(object)
    if (!"sampleNameAggregate" %in% colnames(sampleMetadata)) {
        stop("'sampleNameAggregate' not present in sample metadata")
    }
    # We'll end up replacing `sampleID` and `sampleName` columns with the
    # corresponding `*Aggregate` columns.
    sampleMetadata <-
        sampleMetadata[, c("sampleID",
                           "sampleName",
                           "sampleNameAggregate")] %>%
        mutate(sampleIDAggregate = make.names(
            .data[["sampleNameAggregate"]],
            unique = FALSE)
        ) %>%
        mutate_all(as.factor) %>%
        .[, c("sampleIDAggregate",
              "sampleID",
              "sampleNameAggregate",
              "sampleName")] %>%
        arrange(.data[["sampleIDAggregate"]], .data[["sampleID"]]) %>%
        mutate_all(reorder)

    # Message the new sample IDs
    newIDs <- unique(sampleMetadata[["sampleIDAggregate"]])
    message(paste(
        "New sample IDs:", toString(newIDs)
    ))

    message("Remapping cellular barcodes to aggregate sample IDs")
    cell2sample <- cell2sample(object)
    map <- left_join(
        data.frame(sampleID = cell2sample),
        sampleMetadata,
        by = "sampleID")
    rownames(map) <- names(cell2sample)
    cells <- mapply(
        FUN = gsub,
        x = rownames(map),
        pattern = paste0("^", map[["sampleID"]]),
        replacement = map[["sampleIDAggregate"]]
    )

    # Counts ===================================================================
    message("Aggregating counts")
    counts <- aggregateReplicates(assay(object), cells = cells)
    # Check that the count number of counts matches
    if (!identical(sum(assay(object)), sum(counts))) {
        stop("Aggregated counts sum doens't match the original",
             call. = FALSE)
    }

    # Cellular barcodes list ===================================================
    message("Aggregating raw cellular barcode counts")
    cb <- bcbio(object, "cellularBarcodes")
    if (is.null(cb)) {
        cbAggregateList <- NULL
    } else {
        if (is.list(cb)) {
            # Aggregate and split back out as a list?
            colnames <- c("sampleID", colnames(cb[[1]]))
            cb <- .bindCellularBarcodes(cb)
            cb <- cb[, colnames]
        }
        # Now let's remap the cellular barcode counts for the aggregated samples
        cbAggregateData <- cb %>%
            # Now we need to map the new sampleIDs
            left_join(sampleMetadata[, c("sampleID", "sampleIDAggregate")],
                      by = "sampleID") %>%
            ungroup() %>%
            mutate(sampleID = .data[["sampleIDAggregate"]],
                   sampleIDAggregate = NULL) %>%
            # Here we're grouping per cellular barcode
            group_by(!!!syms(c("sampleID", "cellularBarcode"))) %>%
            # Now sum the counts for each unique barcode
            summarize(nCount = sum(.data[["nCount"]])) %>%
            ungroup() %>%
            group_by(.data[["sampleID"]]) %>%
            arrange(desc(.data[["nCount"]]), .by_group = TRUE)
        # Group and sum the counts
        # Now split this back out into a list to match the original data
        # structure
        cbAggregateList <- lapply(seq_along(newIDs), function(a) {
            cbAggregateData %>%
                ungroup() %>%
                filter(sampleID == newIDs[[a]]) %>%
                mutate(sampleID = NULL)
        })
        names(cbAggregateList) <- as.character(newIDs)
    }

    # Metrics ==================================================================
    annotable <- annotable(object)
    prefilter <- metadata(object)[["prefilter"]]

    # Rather than recalculating, we need to just remap the rownames
    metrics <- calculateMetrics(
        counts,
        annotable = annotable,
        prefilter = prefilter)
    # Add the raw read counts for objects with raw cellular barcodes list
    if (!is.null(cb)) {
        nCount <- cbAggregateData %>%
            as.data.frame() %>%
            ungroup() %>%
            .[, "nCount"]
        rownames(nCount) <- paste(
            cbRemap[["sampleID"]],
            cbRemap[["cellularBarcode"]],
            sep = "_")
        nCount <- nCount[rownames(metrics), , drop = FALSE]
        metrics <- cbind(metrics, nCount)
    }

    # Prefilter very low quality cells, if desired
    if (isTRUE(prefilter)) {
        # Subset the counts matrix to match the metrics
        counts <- counts[, rownames(metrics)]
    }

    # Metadata =================================================================
    message("Updating metadata")
    metadata <- metadata(object)
    metadata[["sampleMetadata"]] <-
        sampleMetadata(object, aggregateReplicates = TRUE)
    metadata[["cell2sample"]] <-
        cell2sample(colnames(counts), samples = newIDs)
    # Slot the named vector used to aggregate the replicates
    metadata[["aggregateReplicates"]] <- cells
    # Update filtered cells. We can use the named `cells` character vector
    # defined above to remap.
    filterCells <- metadata[["filterCells"]]
    if (!is.null(filterCells)) {
        metadata[["filterCells"]] <- cells[filterCells]
    }

    # bcbioSingleCell ==========================================================
    message("Preparing aggregated bcbioSingleCell object")
    se <- SummarizedExperiment(
        assays = list(assay = counts),
        rowData = annotable,
        colData = metrics,
        metadata = metadata)
    bcb <- new("bcbioSingleCell", se)
    bcbio(bcb, "cellularBarcodes") <- cbAggregateList
    validObject(bcb)
    bcb
}



# Methods ====
#' @rdname aggregateReplicates
#' @export
setMethod(
    "aggregateReplicates",
    signature("bcbioSingleCell"),
    .aggregateReplicates)
