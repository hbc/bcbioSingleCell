# FIXME Need to update this method.
# FIXME rownames to column `sampleID`...



#' Aggregate Columns
#'
#' @name aggregateCols
#' @family Data Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @importFrom basejump aggregateCols
#'
#' @inheritParams general
#'
#' @return `SingleCellExperiment`.
#'
#' @examples
#' object <- indrops_small
#' sampleNames(object)
#'
#' # Define groupings factor as`aggregate` column in `colData()`.
#' glimpse(object$aggregate)
#'
#' x <- aggregateCols(object)
#' show(x)
#' sampleNames(x)
NULL



#' @rdname aggregateCols
#' @export
setMethod(
    "aggregateCols",
    signature("SingleCellExperiment"),
    function(object) {
        validObject(object)

        # FIXME Move this sampleData code to SE method?

        # FIXME Move this to SE method.
        # Aggregate the sample data.
        colData <- colData(object)
        if ("sampleNameAggregate" %in% colnames(colData)) {
            stop("Use `aggregate` instead of `sampleNameAggregate`")
        }
        assert_is_subset("aggregate", colnames(colData))

        # FIXME Require valid names in aggregate column?
        # Consider reworking `assertAllAreValidNames()`.
        assert_all_are_true(validNames(
            colData(object)[["aggregate"]]
        ))

        # Consider adding an assert check here to check that interesting
        # groups map to aggregate-level sample columns.
        interestingGroupsAggregate <- setdiff(
            x = interestingGroups(object),
            y = "sampleName"
        )

        # Collapse the sample data ---------------------------------------------
        # This step will replace the `sampleName` column with the `aggregate`
        # column metadata.
        # FIXME Improve the tibble/dplyr code here.
        sampleData <- sampleData(object) %>%
            as("tbl_df") %>%
            select(!!!syms(unique(c(
                "aggregate",
                interestingGroupsAggregate
            )))) %>%
            # FIXME Can take this out if we require aggregate to be valid.
            mutate(rowname = makeNames(
                # Use `unique = TRUE` here instead?
                !!sym("aggregate"), unique = FALSE
            )) %>%
            rename(sampleName = !!sym("aggregate")) %>%
            arrange(!!!syms(c("rowname", "sampleName"))) %>%
            mutate_all(as.factor) %>%
            mutate_all(droplevels) %>%
            as("DataFrame")

        # Message the new sample IDs
        message(paste(
            "New sample names:",
            toString(levels(sampleData[["sampleName"]]))
        ))

        # Remap cellular barcodes ----------------------------------------------
        message("Remapping cellular barcodes to aggregate sample IDs")
        cell2sample <- cell2sample(object)
        map <- tibble(
            cell = names(cell2sample),
            sample = cell2sample,
            aggregate = colData(object)[["aggregate"]]
        )
        groupings <- mapply(
            FUN = gsub,
            x = map[["cell"]],
            pattern = paste0("^", map[["sample"]]),
            replacement = map[["aggregate"]],
            SIMPLIFY = TRUE,
            USE.NAMES = TRUE
        ) %>%
            as.factor()
        # Reslot the `aggregate` column using these groupings.
        assert_are_identical(names(groupings), colnames(object))
        colData(object)[["aggregate"]] <- groupings

        # Coerce to RangedSummarizedExperiment and aggregate.
        agg <- aggregateCols(as(object, "RangedSummarizedExperiment"))

        # Column data ----------------------------------------------------------
        # FIXME Rework this?

        # Cell to sample mappings
        cell2sample <- .mapCellsToSamples(
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
        cell2sample <- cell2sample[colnames(counts)]

        # Metadata -------------------------------------------------------------
        metadata <- list(
            aggregateCols = groupings,
            cell2sample = cell2sample,
            interestingGroups = interestingGroups(object),
            sampleData = sampleData
        )

        # Return ---------------------------------------------------------------
        .new.SingleCellExperiment(
            assays = list(counts = counts),
            rowRanges = rowRanges(object),
            colData = colData,
            metadata = metadata,
            spikeNames = spikeNames(object)
        )
    }
)
