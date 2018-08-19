#' Barcode Ranks per Sample
#'
#' @name barcodeRanksPerSample
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inherit DropletUtils::barcodeRanks
#'
#' @inheritParams general
#' @param ... Passthrough arguments to [DropletUtils::barcodeRanks()].
#'
#' @seealso [DropletUtils::barcodeRanks()].
#'
#' @examples
#' x <- barcodeRanksPerSample(indrops_small)
#' names(x)
NULL



#' @rdname barcodeRanksPerSample
#' @export
setMethod(
    "barcodeRanksPerSample",
    signature("SingleCellExperiment"),
    function(object, ...) {
        # DropletUtils requires R 3.5, which isn't on conda.
        # Don't import the package as a dependency for the time being.
        requireNamespace("DropletUtils")
        counts <- counts(object)
        cell2sample <- cell2sample(object)
        samples <- levels(cell2sample)
        perSampleCounts <- lapply(samples, function(sample) {
            cells <- names(cell2sample)[which(cell2sample == sample)]
            counts[, cells]
        })
        ranks <- lapply(perSampleCounts, function(object) {
            DropletUtils::barcodeRanks(as.matrix(object), ...)
        })
        names(ranks) <- samples
        ranks
    }
)
