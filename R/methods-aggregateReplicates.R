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

        # Require that the dataset is filtered
        if (!length(metadata(object)[["filterCells"]])) {
            stop(paste(
                "`aggregateReplicates()` is only supported for filtered",
                "bcbioSingleCell objects.",
                "Run `filterCells()` prior to replicate aggregation."
            ))
        }

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
            mutate_all(droplevels)

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
        counts <- aggregateReplicates(counts(object), groupings = groupings)
        # Check that the count number of counts matches
        if (!identical(sum(assay(object)), sum(counts))) {
            stop("Aggregated counts sum isn't identical to original")
        }

        # Row data =============================================================
        rowData <- rowData(object)
        rownames(rowData) <- rownames(object)

        # Column data ==========================================================
        # Always prefilter, removing cells with no UMIs or genes
        metrics <- metrics(counts, rowData = rowData, prefilter = TRUE)

        # Cell to sample mappings
        cell2sample <- mapCellsToSamples(
            cells = rownames(metrics),
            samples = rownames(sampleData)
        )

        sampleData[["sampleID"]] <- rownames(sampleData)
        colData <- as(metrics, "DataFrame")
        colData[["cellID"]] <- rownames(colData)
        colData[["sampleID"]] <- cell2sample
        colData <- merge(
            x = colData,
            y = sampleData,
            by = "sampleID",
            all.x = TRUE
        )
        rownames(colData) <- colData[["cellID"]]
        colData[["cellID"]] <- NULL
        sampleData[["sampleID"]] <- NULL

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

        # filterCells
        filterCells <- metadata[["filterCells"]]
        assert_is_character(filterCells)
        metadata[["filterCells"]] <- groupings[filterCells]

        # aggregateReplicates
        metadata[["aggregateReplicates"]] <- groupings

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
