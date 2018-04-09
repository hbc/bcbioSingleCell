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
        objects <- subsetPerSample(object, assignAndSave = FALSE)
        lapply(objects, barcodeRanks)
    }
)
