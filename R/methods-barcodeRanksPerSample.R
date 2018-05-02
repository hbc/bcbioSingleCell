#' Barcode Ranks per Sample
#'
#' @name barcodeRanksPerSample
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inherit barcodeRanks
#' @param ... Passthrough arguments to [barcodeRanks()].
#'
#' @examples
#' # bcbioSingleCell ====
#' x <- barcodeRanksPerSample(indrops_small)
#' names(x)
#'
#' # seurat ====
#' x <- barcodeRanksPerSample(seurat_small)
#' names(x)
NULL



# Methods ======================================================================
#' @rdname barcodeRanksPerSample
#' @export
setMethod(
    "barcodeRanksPerSample",
    signature("SingleCellExperiment"),
    function(object, ...) {
        counts <- counts(object)
        cell2sample <- cell2sample(object)
        samples <- levels(cell2sample)
        perSampleCounts <- lapply(samples, function(sample) {
            cells <- names(cell2sample)[which(cell2sample == sample)]
            counts[, cells]
        })
        ranks <- lapply(perSampleCounts, function(object) {
            barcodeRanks(object, ...)
        })
        names(ranks) <- samples
        ranks
    }
)



#' @rdname barcodeRanksPerSample
#' @export
setMethod(
    "barcodeRanksPerSample",
    signature("seurat"),
    getMethod("barcodeRanksPerSample", "SingleCellExperiment")
)
