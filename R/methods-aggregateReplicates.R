#' Aggregate Replicates
#'
#' @name aggregateReplicates
#' @family Data Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom basejump aggregateReplicates
#'
#' @inheritParams general
#'
#' @return `bcbioSingleCell`.
#'
#' @examples
#' # bcbioSingleCell ====
#' x <- aggregateReplicates(indrops_small)
#' sampleData(x)
NULL



# Methods ======================================================================
#' @rdname aggregateReplicates
#' @export
setMethod(
    "aggregateReplicates",
    signature("bcbioSingleCell"),
    function(object) {
        validObject(object)
        metadata <- metadata(object)
        sampleData <- sampleData(object, clean = FALSE, return = "data.frame")
        if ("sampleNameAggregate" %in% colnames(sampleData)) {
            warning("Use `aggregate` instead of `sampleNameAggregate`")
            sampleData[["aggregate"]] <- sampleData[["sampleNameAggregate"]]
        }
        assert_is_subset("aggregate", colnames(sampleData))

        # This step will replace the `sampleName` column with the `aggregate`
        # column metadata.
        map <- sampleData %>%
            rownames_to_column("sampleID") %>%
            select(!!!syms(c("sampleID", "aggregate"))) %>%
            mutate(sampleIDAggregate = makeNames(
                !!sym("aggregate"), unique = FALSE
            )) %>%
            select(-!!sym("aggregate")) %>%
            arrange(!!!syms(c("sampleID", "sampleIDAggregate"))) %>%
            mutate_all(as.factor) %>%
            mutate_all(reorder)

        # Message the new sample IDs
        newIDs <- unique(map[["sampleIDAggregate"]])
        message(paste("New sample IDs:", toString(newIDs)))

        message("Remapping cellular barcodes to aggregate sample IDs")
        cell2sample <- cell2sample(object)
        remap <- tibble(
            "cellID" = names(cell2sample),
            "sampleID" = cell2sample
        ) %>%
            left_join(map, by = "sampleID")

        groupings <- mapply(
            FUN = gsub,
            x = remap[["cellID"]],
            pattern = paste0("^", remap[["sampleID"]]),
            replacement = remap[["sampleIDAggregate"]]
        ) %>%
            as.factor()

        # Assays ===============================================================
        message("Aggregating counts")
        counts <- aggregateReplicates(assay(object), groupings = groupings)
        # Check that the count number of counts matches
        if (!identical(sum(assay(object)), sum(counts))) {
            stop("Aggregated counts sum isn't identical to original")
        }

        # Row data =============================================================
        rowData <- rowData(object)
        rownames(rowData) <- rownames(object)

        # Column data ==========================================================
        # Always prefilter, removing cells with no UMIs or genes
        colData <- metrics(counts, rowData = rowData, prefilter = TRUE)

        # Subset the counts to match the prefiltered metrics
        counts <- counts[, rownames(colData), drop = FALSE]

        # Metadata =============================================================
        message("Updating metadata")

        # sampleData
        expected <- length(levels(sampleData[["aggregate"]]))
        sampleData <- sampleData %>%
            mutate(sampleName = !!sym("aggregate")) %>%
            select(!!sym("sampleName")) %>%
            mutate_all(as.factor) %>%
            unique()
        if (!identical(nrow(sampleData), expected)) {
            stop("Failed to aggregate sample metadata uniquely")
        }
        rownames(sampleData) <- makeNames(sampleData[["sampleName"]])
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
            cb <- cb[, colnames, drop = FALSE]
            # Now let's remap the barcode reads for the aggregated samples
            cbAggregateData <- cb %>%
                # Now we need to map the new sampleIDs
                left_join(
                    map[, c("sampleID", "sampleIDAggregate")],
                    by = "sampleID"
                ) %>%
                ungroup() %>%
                mutate(
                    sampleID = !!sym("sampleIDAggregate"),
                    sampleIDAggregate = NULL
                ) %>%
                # Here we're grouping per cellular barcode
                group_by(!!!syms(c("sampleID", "cellularBarcode"))) %>%
                # Now sum the counts for each unique barcode.
                # This step is CPU intensive when there's many samples and we
                # may want to add a progress bar here in a future update.
                summarize(nCount = sum(!!sym("nCount"))) %>%
                ungroup() %>%
                group_by(!!sym("sampleID")) %>%
                arrange(desc(!!sym("nCount")), .by_group = TRUE)
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

        # Return ===============================================================
        .new.bcbioSingleCell(
            assays = list("counts" = counts),
            rowRanges = rowRanges(object),
            colData = colData,
            metadata = metadata,
            spikeNames = spikeNames(object)
        )
    }
)
