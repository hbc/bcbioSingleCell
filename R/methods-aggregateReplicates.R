#' Aggregate Replicates
#'
#' @name aggregateReplicates
#' @family Data Management Utilities
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom basejump aggregateReplicates
#'
#' @inheritParams general
#'
#' @examples
#' # bcbioSingleCell ====
#' x <- aggregateReplicates(bcb_small)
#' sampleData(x)
NULL



# Constructors =================================================================
#' @importFrom basejump makeNames
#' @importFrom dplyr arrange group_by mutate mutate_all ungroup
#' @importFrom parallel mclapply
#' @importFrom rlang !!! syms
#' @importFrom stats reorder
#' @importFrom stringr str_match
#' @importFrom tibble column_to_rownames rownames_to_column
.aggregateReplicates <- function(object) {
    validObject(object)
    metadata <- metadata(object)

    sampleData <- sampleData(object)
    assert_is_subset("sampleNameAggregate", colnames(sampleData))
    # We'll end up replacing `sampleID` and `sampleName` columns with the
    # corresponding `*Aggregate` columns.
    map <- sampleData %>%
        .[, c("sampleID", "sampleName", "sampleNameAggregate")] %>%
        mutate(sampleIDAggregate = makeNames(
            .data[["sampleNameAggregate"]],
            unique = FALSE
        )) %>%
        mutate_all(as.factor) %>%
        .[, c(
            "sampleIDAggregate", "sampleID",
            "sampleNameAggregate",
            "sampleName"
        )] %>%
        arrange(.data[["sampleIDAggregate"]], .data[["sampleID"]]) %>%
        mutate_all(reorder)

    # Message the new sample IDs
    newIDs <- unique(map[["sampleIDAggregate"]])
    inform(paste(
        "New sample IDs:", toString(newIDs)
    ))

    inform("Remapping cellular barcodes to aggregate sample IDs")
    cell2sample <- cell2sample(object)
    sampleID <- data.frame(sampleID = cell2sample)
    remap <- left_join(
        x = sampleID,
        y = map,
        by = "sampleID"
    )
    rownames(remap) <- names(cell2sample)
    groupings <- mapply(
        FUN = gsub,
        x = rownames(remap),
        pattern = paste0("^", remap[["sampleID"]]),
        replacement = remap[["sampleIDAggregate"]]
    ) %>%
        as.factor()

    # Assays ===================================================================
    inform("Aggregating counts")
    counts <- aggregateReplicates(assay(object), groupings = groupings)
    # Check that the count number of counts matches
    if (!identical(sum(assay(object)), sum(counts))) {
        abort("Aggregated counts sum doens't match the original")
    }

    # Column data ==============================================================
    rowData <- rowData(object)
    prefilter <- metadata[["prefilter"]]

    # Rather than recalculating, we need to just remap the rownames
    colData <- calculateMetrics(
        counts,
        rowData = rowData,
        prefilter = prefilter
    )

    # Prefilter very low quality cells, if desired
    if (isTRUE(prefilter)) {
        # Subset the counts matrix to match the colData
        counts <- counts[, rownames(colData)]
    }

    # Metadata =================================================================
    inform("Updating metadata")

    # sampleData
    expected <- length(unique(sampleData[["sampleNameAggregate"]]))
    sampleData <- sampleData %>%
        mutate(
            sampleName = .data[["sampleNameAggregate"]],
            description = .data[["sampleName"]],
            sampleID = makeNames(
                .data[["sampleName"]], unique = FALSE
            )
        ) %>%
        # TODO Improve detection and handling of unique columns here
        .[, metadataPriorityCols] %>%
        unique()
    if (!identical(nrow(sampleData), expected)) {
        abort("Failed to aggregate sample metadata uniquely")
    }
    rownames(sampleData) <- sampleData[["sampleID"]]
    metadata[["sampleData"]] <- sampleData

    # cell2sample
    metadata[["cell2sample"]] <- mapCellsToSamples(
        cells = colnames(counts),
        samples = newIDs
    )

    # aggregateReplicates
    metadata[["aggregateReplicates"]] <- groupings

    # cellularBarcodes, if set
    cb <- metadata[["cellularBarcodes"]]
    if (is.list(cb)) {
        # Aggregate and split back out as a list?
        colnames <- c("sampleID", colnames(cb[[1L]]))
        cb <- .bindCellularBarcodes(cb)
        cb <- cb[, colnames]
        # Now let's remap the cellular barcode counts for the aggregated samples
        cbAggregateData <- cb %>%
            # Now we need to map the new sampleIDs
            left_join(
                map[, c("sampleID", "sampleIDAggregate")],
                by = "sampleID"
            ) %>%
            ungroup() %>%
            mutate(
                sampleID = .data[["sampleIDAggregate"]],
                sampleIDAggregate = NULL
            ) %>%
            # Here we're grouping per cellular barcode
            group_by(!!!syms(c("sampleID", "cellularBarcode"))) %>%
            # Now sum the counts for each unique barcode.
            # This step is CPU intensive when there's many samples and we
            # may want to add a progress bar here in a future update.
            summarize(nCount = sum(.data[["nCount"]])) %>%
            ungroup() %>%
            group_by(.data[["sampleID"]]) %>%
            arrange(dplyr::desc(.data[["nCount"]]), .by_group = TRUE)
        # Group and sum the counts
        # Now split this back out into a list to match the original data
        # structure
        cbAggregateList <- lapply(seq_along(newIDs), function(a) {
            cbAggregateData %>%
                ungroup() %>%
                .[.[["sampleID"]] == newIDs[[a]], , drop = FALSE] %>%
                mutate(sampleID = NULL)
        })
        names(cbAggregateList) <- as.character(newIDs)
        metadata[["cellularBarcodes"]] <- cbAggregateList
    }

    # filterCells, if set
    filterCells <- metadata[["filterCells"]]
    if (!is.null(filterCells)) {
        metadata[["filterCells"]] <- groupings[filterCells]
    }

    # Return ===================================================================
    .new.bcbioSingleCell(
        assays = list(assay = counts),
        rowRanges = rowRanges(object),
        colData = colData,
        metadata = metadata,
        isSpike = isSpike(object)
    )
}



# Methods ======================================================================
#' @rdname aggregateReplicates
#' @export
setMethod(
    "aggregateReplicates",
    signature("bcbioSingleCell"),
    .aggregateReplicates
)
