#' Barcode Ranks per Sample
#'
#' @name barcodeRanksPerSample
#' @author Michael Steinbaugh
#'
#' @inheritParams barcodeRanks
#'
#' @examples
#' # bcbioSingleCell ====
#' x <- barcodeRanksPerSample(bcb_small)
#' names(x)
NULL



# Methods ======================================================================
#' @rdname barcodeRanksPerSample
#' @export
setMethod(
    "barcodeRanksPerSample",
    signature("SingleCellExperiment"),
    function(object) {
        counts <- counts(object)
        cell2sample <- cell2sample(object)
        samples <- levels(cell2sample)
        perSampleCounts <- lapply(samples, function(sample) {
            cells <- names(cell2sample)[which(cell2sample == sample)]
            counts[, cells]
        })
        ranks <- lapply(perSampleCounts, barcodeRanks)
        names(ranks) <- samples
        ranks
    }
)
