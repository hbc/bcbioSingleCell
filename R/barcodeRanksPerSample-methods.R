#' @name barcodeRanksPerSample
#' @author Michael Steinbaugh
#' @inherit bioverbs::barcodeRanksPerSample
#' @inherit DropletUtils::barcodeRanks
#' @inheritParams acidplots::params
#' @inheritParams basejump::params
#' @param ... Additional arguments.
#' @seealso [DropletUtils::barcodeRanks()].
#' @examples
#' data(indrops)
#' if (packageVersion("DropletUtils") >= "1.4") {
#'     x <- barcodeRanksPerSample(indrops)
#'     names(x)
#' }
NULL



#' @rdname barcodeRanksPerSample
#' @name barcodeRanksPerSample
#' @importFrom bioverbs barcodeRanksPerSample
#' @usage barcodeRanksPerSample(object, ...)
#' @export
NULL



barcodeRanksPerSample.bcbioSingleCell <-  # nolint
    function(object) {
        assert(packageVersion("DropletUtils") >= "1.4")
        which <- sys.parent()

        counts <- counts(object)
        cell2sample <- cell2sample(object)
        samples <- levels(cell2sample)

        ## Subset the counts per sample into a list.
        countsPerSample <- lapply(samples, function(sample) {
            cells <- names(cell2sample)[which(cell2sample == sample)]
            counts[, cells, drop = FALSE]
        })

        ## Calculate the ranks per sample.
        ranks <- lapply(
            X = countsPerSample,
            FUN = function(counts) {
                x <- do.call(
                    what = barcodeRanks,
                    args = matchArgsToDoCall(
                        args = list(m = as.matrix(counts)),
                        removeFormals = "object",
                        which = which
                    )
                )

                ## Check DropletUtils return.
                assert(
                    is(x, "DataFrame"),
                    identical(
                        x = colnames(x),
                        y = c("rank", "total", "fitted")
                    ),
                    isSubset(
                        x = names(metadata(x)),
                        y = c("knee", "inflection")
                    )
                )

                x
            }
        )

        names(ranks) <- samples
        ranks
    }

f1 <- formals(barcodeRanksPerSample.bcbioSingleCell)
f2 <- formals(barcodeRanks)
f2 <- f2[setdiff(names(f2), c(names(f1), "m", "..."))]
f <- c(f1, f2)
formals(barcodeRanksPerSample.bcbioSingleCell) <- f



#' @rdname barcodeRanksPerSample
#' @export
setMethod(
    f = "barcodeRanksPerSample",
    signature = signature("bcbioSingleCell"),
    definition = barcodeRanksPerSample.bcbioSingleCell
)
