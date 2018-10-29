#' Barcode Ranks per Sample
#'
#' @name barcodeRanksPerSample
#' @author Michael Steinbaugh
#'
#' @inherit DropletUtils::barcodeRanks
#'
#' @inheritParams basejump.globals::params
#' @param ... Additional arguments.
#'
#' @seealso [DropletUtils::barcodeRanks()].
#'
#' @examples
#' data(indrops_small)
#' x <- barcodeRanksPerSample(indrops_small)
#' names(x)
NULL



barcodeRanksPerSample.SingleCellExperiment <-  # nolint
    function(object) {
        which <- sys.parent()

        counts <- counts(object)
        cell2sample <- cell2sample(object)
        samples <- levels(cell2sample)

        # Subset the counts per sample into a list.
        countsPerSample <- lapply(samples, function(sample) {
            cells <- names(cell2sample)[which(cell2sample == sample)]
            counts[, cells, drop = FALSE]
        })

        # Calculate the ranks per sample.
        ranks <- lapply(
            X = countsPerSample,
            FUN = function(counts) {
                do.call(
                    what = barcodeRanks,
                    args = matchArgsToDoCall(
                        args = list(m = as.matrix(counts)),
                        removeFormals = "object",
                        which = which
                    )
                )
            }
        )

        names(ranks) <- samples
        ranks
    }
f1 <- formals(barcodeRanksPerSample.SingleCellExperiment)
f2 <- formals(barcodeRanks)
f2 <- f2[setdiff(names(f2), c(names(f1), "m", "..."))]
f <- c(f1, f2)
formals(barcodeRanksPerSample.SingleCellExperiment) <- f



#' @rdname barcodeRanksPerSample
#' @export
setMethod(
    f = "barcodeRanksPerSample",
    signature = signature("SingleCellExperiment"),
    definition = barcodeRanksPerSample.SingleCellExperiment
)
