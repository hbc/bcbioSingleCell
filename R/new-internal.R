.new.SingleCellExperiment <- function(  # nolint
    assays,
    rowRanges,
    colData,
    metadata,
    transgeneNames = NULL,
    spikeNames = NULL
) {
    # Prepare RangedSummarizedExperiment first.
    # Supports automatic resizing of rowRanges and helps slot FASTA spike-ins.
    rse <- prepareSummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )
    # Then coerce to SingleCellExperiment
    sce <- SingleCellExperiment(
        assays = assays(rse),
        rowRanges = rowRanges(rse),
        colData = colData(rse),
        metadata = metadata(rse)
    )

    # Optionally, use `isSpike` internally to define the `spikeNames`
    if (is.character(spikeNames)) {
        for (i in seq_along(spikeNames)) {
            isSpike(sce, spikeNames[[i]]) <- spikeNames[[i]]
        }
    }

    sce
}



.new.bcbioSingleCell <-  # nolint
    function(...) {
        sce <- .new.SingleCellExperiment(...)
        new("bcbioSingleCell", sce)
    }
