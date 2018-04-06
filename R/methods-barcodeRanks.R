# TODO Import DropletUtils when Bioc 3.7 is released

#' Barcode Ranks
#'
#' Return statistics to construct a barcode-rank plot. Also calculates the knee
#' and inflection points for further use.
#'
#' @name barcodeRanks
#' @author Aaron Lun, Michael Steinbaugh
#'
#' @seealso Modified version of `DropletUtils::barcodeRanks()` that works on
#'   a `SingleCellExperiment` object instead of a matrix.
#'
#' @return `list`.
#' @examples
#' # SingleCellExperiment ====
#' x <- barcodeRanks(bcb_small)
#' names(x)
#' x[["knee"]]
#' x[["inflection"]]
NULL



# Methods ======================================================================
#' @rdname barcodeRanks
#' @export
setMethod(
    "barcodeRanks",
    signature("SingleCellExperiment"),
    function(
        object,
        lower = 100L,
        fitBounds = NULL,
        df = 20L
    ) {
        m <- counts(object)
        totals <- Matrix::colSums(m)
        o <- order(totals, decreasing = TRUE)

        # Run length encoding
        rle <- rle(totals[o])

        # Get mid-rank of each run
        runRank <- cumsum(rle[["lengths"]]) - (rle[["lengths"]] - 1L) / 2L
        runTotals <- rle[["values"]]

        keep <- runTotals > lower
        if (sum(keep) < 3L) {
            abort("Insufficient totals for computing knee/inflection points")
        }
        y <- log10(runTotals[keep])
        x <- log10(runRank[keep])

        # Numerical differentiation to identify bounds for spline fitting.
        # The upper/lower bounds are defined at the plateau and inflection.
        d1n <- diff(y) / diff(x)
        rightEdge <- which.min(d1n)
        leftEdge <- which.max(d1n[seq_len(rightEdge)])

        # We restrict to this region, thereby simplifying the shape of the
        # curve. This allows us to get a decent fit with low df for stable
        # differentiation.
        if (is.null(fitBounds)) {
            newKeep <- leftEdge:rightEdge
        } else {
            newKeep <- y > log10(fitBounds[[1L]]) & y < log10(fitBounds[[2L]])
        }

        # Smoothing to avoid error multiplication upon differentiation.
        # Minimizing the signed curvature and returning the total for the knee
        # point.
        fit <- smooth.spline(
            x = x[newKeep],
            y = y[newKeep],
            df = df
        )

        # Calculate derivatives
        d1 <- predict(fit, deriv = 1L)[["y"]]
        d2 <- predict(fit, deriv = 2L)[["y"]]
        curvature <- d2 / (1L + d1 ^ 2L) ^ 1.5

        # Knee point
        knee <- 10L ^ (y[which.min(curvature)])

        # Taking the right edge to get the total for the inflection point. We
        # use the numerical derivative as the spline is optimized for the knee.
        inflection <- 10L ^ (y[rightEdge])

        # Returning a whole stack of useful stats
        fittedVals <- rep(NA_real_, length(keep))
        fittedVals[keep][newKeep] <- 10L ^ fitted(fit)

        list(
            rank = .reorder(runRank, rle$lengths, o),
            total = .reorder(runTotals, rle$lengths, o),
            fitted = .reorder(fittedVals, rle$lengths, o),
            knee = knee,
            inflection = inflection
        )
    }
)

.reorder <- function(vals, lens, o) {
    out <- rep(vals, lens)
    out[o] <- out
    out
}
