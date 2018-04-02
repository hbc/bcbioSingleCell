#' Calculate Knee Plot Inflection Point
#'
#' @name inflectionPoint
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#' @param maxCells The number of cells to consider as an upper bound and to
#'   compute the inflection point.
#'
#' @seealso Modified version of `dropbead::estimateCellNumber()`.
#'
#' @examples
#' inflectionPoint(bcb_small)
NULL



# Constructors =================================================================
.inflectionPoint <- function(
    object,
    maxCells = 100000L
) {
    validObject(object)

    # Get the read counts
    metrics <- metrics(object)
    countsCol <- c("nCount", "nUMI")
    assert_are_intersecting_sets(countsCol, colnames(metrics))
    countsCol <- intersect(countsCol, colnames(metrics))[[1L]]
    counts <- metrics[[countsCol]]

    # Apply the upper boundary cutoff defined by `maxCells`
    cutoff <- min(length(counts), maxCells)
    counts <- counts[seq_len(cutoff)]

    # Apply trigonometric calculations
    xdata <- seq_along(counts) / length(counts)
    cs <- cumsum(counts / sum(counts))
    m <- (cs[[cutoff]] - cs[[1L]]) / (cutoff - xdata[[1L]])
    dists <- vapply(
        X = seq_len(cutoff),
        FUN = function(cell) {
            sin(atan(m)) * abs(cs[[cell]] - xdata[[cell]])
        },
        FUN.VALUE = numeric(1L)
    )

    # Return the count cutoff threshold (not the cell/row number)
    counts[which.max(dists)]
}



# Methods ======================================================================
#' @rdname inflectionPoint
#' @export
setMethod(
    "inflectionPoint",
    signature("bcbioSingleCell"),
    .inflectionPoint
)
