#' Barcode Ranks per Sample
#'
#' @name barcodeRanksPerSample
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#' @export
#'
#' @inherit DropletUtils::barcodeRanks
#'
#' @inheritParams general
#' @param ... Additional arguments.
#'
#' @seealso [DropletUtils::barcodeRanks()].
#'
#' @examples
#' x <- barcodeRanksPerSample(indrops_small)
#' names(x)
NULL



.barcodeRanksPerSample <- function(object) {
    # Match the call, which we need to pass to barcodeRanks below.
    call <- matchCall()
    fun <- sys.function(which = sys.parent())

    counts <- counts(object)
    cell2sample <- cell2sample(object)
    samples <- levels(cell2sample)

    # Subset the counts per sample into a list.
    countsPerSample <- lapply(samples, function(sample) {
        cells <- names(cell2sample)[which(cell2sample == sample)]
        counts[, cells, drop = FALSE]
    })

    # Calculate the ranks per sample.
    ranks <- lapply(countsPerSample, function(object) {
        args <- setArgsToDoCall(
            args = list(m = as.matrix(object)),
            removeArgs = "object",
            call = call,
            fun = fun
        )
        do.call(what = barcodeRanks, args = args)
    })

    names(ranks) <- samples
    ranks
}

# Assign the formals.
f1 <- formals(.barcodeRanksPerSample)
f2 <- formals(barcodeRanks)
f2 <- f2[setdiff(names(f2), c(names(f1), "m", "..."))]
f <- c(f1, f2)
formals(.barcodeRanksPerSample) <- f



#' @rdname barcodeRanksPerSample
#' @export
setMethod(
    f = "barcodeRanksPerSample",
    signature = signature("SingleCellExperiment"),
    definition = .barcodeRanksPerSample
)
