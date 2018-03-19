#' @importFrom bcbioBase prepareSummarizedExperiment
#' @importFrom SingleCellExperiment SingleCellExperiment
.newBcbioSingleCell <- function(
    assays,
    rowRanges,
    colData,
    metadata,
    isSpike
) {
    metadata <- Filter(Negate(is.null), metadata)

    # Prepare RangedSummarizedExperiment, with automatic resizing of rowRanges
    # and support for FASTA spike-ins
    rse <- prepareSummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata,
        isSpike = isSpike
    )

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

    new("bcbioSingleCell", sce)
}
