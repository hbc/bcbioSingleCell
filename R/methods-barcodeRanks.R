#' Barcode Ranks
#'
#' Return statistics to construct a barcode-rank plot. Also calculates the knee
#' and inflection points for further use.
#'
#' @name barcodeRanks
#' @family Quality Control Functions
#' @author Aaron Lun, Michael Steinbaugh
#'
#' @inheritParams general
#' @param lower `integer` scalar specifying the lower bound on the total UMI
#'   count, at or below which all barcodes are assumed to correspond to empty
#'   droplets.
#' @param fitBounds `numeric` vector of length 2, specifying the lower and upper
#'   bounds on the total UMI count for spline fitting.
#' @param df `integer` scalar specifying the number of degress of freedom to
#'   use for spline fitting with [stats::smooth.spline()].
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



# Constructors =================================================================
.reorder <- function(vals, lens, o) {
    assert_has_names(vals)
    out <- rep(vals, lens)
    out[o] <- out
    out
}



# Methods ======================================================================
#' @rdname barcodeRanks
#' @export
setMethod(
    "barcodeRanks",
    signature("matrix"),
    function(
        object,
        lower = 100L,
        fitBounds = NULL,
        df = 20L
    ) {
        assertIsAnImplicitInteger(lower)
        assert_is_any_of(fitBounds, c("numeric", "NULL"))
        if (is.numeric(fitBounds)) {
            assert_is_of_length(fitBounds, 2L)
        }
        assertIsAnImplicitInteger(df)

        # Calculate the total number of UMIs per cell (nUMI)
        totals <- colSums(object)

        # Run length encoding
        o <- order(totals, decreasing = TRUE)
        rle <- rle(totals[o])

        # Get mid-rank of each run
        runTotals <- rle[["values"]]
        runRank <- cumsum(rle[["lengths"]]) - (rle[["lengths"]] - 1L) / 2L

        # Apply lower limit
        keep <- runTotals > lower
        if (sum(keep) < 3L) {
            stop("Insufficient points for computing knee/inflection")
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
        fitNames <- names(y[newKeep])

        # Fitted values
        fittedVals <- rep(NA_real_, length(keep))
        fittedVals[keep][newKeep] <- 10L ^ fitted(fit)
        names(fittedVals) <- names(keep)

        # Calculate derivatives
        d1 <- predict(fit, deriv = 1L)[["y"]]
        names(d1) <- fitNames
        d2 <- predict(fit, deriv = 2L)[["y"]]
        names(d2) <- fitNames
        curvature <- d2 / (1L + d1 ^ 2L) ^ 1.5
        names(curvature) <- fitNames

        # Knee point
        knee <- 10L ^ (y[which.min(curvature)])

        # Take the right edge to get the total for the inflection point.
        # Use the numerical derivative as the spline is optimized for the knee.
        inflection <- 10L ^ (y[rightEdge])

        list(
            "sum" = .reorder(runTotals, rle[["lengths"]], o),
            "midrank" = .reorder(runRank, rle[["lengths"]], o),
            "fitted" = .reorder(fittedVals, rle[["lengths"]], o),
            "fit" = fit,
            "d1" = d1,
            "d2" = d2,
            "curvature" = curvature,
            "knee" = knee,
            "inflection" = inflection
        )
    }
)



#' @rdname barcodeRanks
#' @export
setMethod(
    "barcodeRanks",
    signature("dgCMatrix"),
    getMethod("barcodeRanks", "matrix")
)



#' @rdname barcodeRanks
#' @export
setMethod(
    "barcodeRanks",
    signature("dgTMatrix"),
    getMethod("barcodeRanks", "matrix")
)



#' @rdname barcodeRanks
#' @export
setMethod(
    "barcodeRanks",
    signature("SingleCellExperiment"),
    function(object, ...) {
        validObject(object)
        counts <- counts(object)
        barcodeRanks(counts, ...)
    }
)
