#' Aggregate Columns or Rows
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



# aggregateCols ================================================================
#' @rdname aggregateCols
#' @export
setMethod(
    "aggregateCols",
    signature("SingleCellExperiment"),
    function(object) {
        validObject(object)

        # Remap cellular barcode groupings -------------------------------------
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
        # Add the new cell mappings to the tibble.
        map[["newCell"]] <- groupings
        # Reslot the `aggregate` column using these groupings.
        assert_are_identical(names(groupings), colnames(object))
        colData(object)[["aggregate"]] <- groupings

        # Generate SingleCellExperiment ----------------------------------------
        # Using `SummarizedExperiment` method here.
        rse <- aggregateCols(as(object, "RangedSummarizedExperiment"))
        assert_is_all_of(rse, "RangedSummarizedExperiment")

        # Update the cell-to-sample mappings.
        cell2sample <- .mapCellsToSamples(
            cells = as.character(map[["newCell"]]),
            samples = as.character(map[["aggregate"]])
        )

        # Update the sample data.
        colData(rse)[["sampleID"]] <- cell2sample
        colData(rse)[["sampleName"]] <- cell2sample

        # Now ready to generate aggregated SCE.
        sce <- .new.SingleCellExperiment(
            assays = assays(rse),
            rowRanges = rowRanges(object),
            colData = colData(rse),
            metadata = list(
                aggregateCols = groupings,
                cell2sample = cell2sample,
                interestingGroups = interestingGroups(object)
            ),
            spikeNames = spikeNames(object)
        )

        # Recalculate the metrics.
        sce <- metrics(sce, recalculate = TRUE)

        sce
    }
)
