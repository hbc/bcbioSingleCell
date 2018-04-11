.new.SingleCellExperiment <- function(  # nolint
    assays,
    rowRanges,
    colData,
    metadata,
    transgeneNames,
    spikeNames
) {
    # Prepare RangedSummarizedExperiment. Supports automatic resizing of
    # rowRanges and helps slot FASTA spike-ins.
    rse <- prepareSummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        transgeneNames = transgeneNames,
        spikeNames = spikeNames
    )

    # Coerce to SingleCellExperiment
    sce <- SingleCellExperiment(
        assays = assays(rse),
        rowRanges = rowRanges(rse),
        colData = colData(rse),
        metadata = metadata(rse)
    )

    # Define spikeNames for spike-in sequences
    if (is.character(isSpike)) {
        for (i in seq_along(isSpike)) {
            isSpike(sce, isSpike[[i]]) <- isSpike[[i]]
        }
    }

    sce
}



.new.bcbioSingleCell <-  # nolint
    function(...) {
        sce <- .new.SingleCellExperiment(...)
        new("bcbioSingleCell", sce)
    }
