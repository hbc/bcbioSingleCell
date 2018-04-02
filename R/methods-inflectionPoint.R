#' Calculate the Cellular Barcode Read Depth Inflection Point
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
#' # bcbioSingleCell ====
#' inflectionPoint(bcb_small)
#' inflectionPoint(cellranger_small)
#'
#' # seurat ====
#' inflectionPoint(seurat_small)
#' inflectionPoint(pbmc_small)
NULL



# Constructors =================================================================
.inflectionPoint <- function(
    object,
    maxCells = 100000L
) {
    validObject(object)

    # Get the read counts. Use the raw reads by default. Fall back to plotting
    # UMI count distributions for cellranger and seurat.
    metrics <- metrics(object)
    countsCol <- c("nCount", "nUMI")
    assert_are_intersecting_sets(countsCol, colnames(metrics))
    countsCol <- intersect(countsCol, colnames(metrics))[[1L]]
    counts <- metrics[[countsCol]]

    # Define the upper boundary cutoff, gated by `maxCells` if necessary
    maxCells <- min(length(counts), maxCells)
    counts <- counts[seq_len(maxCells)]

    # Trigonometric calculations
    xdata <- seq_along(counts) / length(counts)
    cs <- cumsum(counts / sum(counts))
    m <- (cs[[maxCells]] - cs[[1L]]) / (maxCells - xdata[[1L]])
    dists <- vapply(
        X = seq_len(maxCells),
        FUN = function(cell) {
            sin(atan(m)) * abs(cs[[cell]] - xdata[[cell]])
        },
        FUN.VALUE = numeric(1L)
    )

    # Return the inflection point as the expression value
    inflection <- counts[which.max(dists)]
    names(inflection) <- countsCol
    inflection
}



# Methods ======================================================================
#' @rdname inflectionPoint
#' @export
setMethod(
    "inflectionPoint",
    signature("bcbioSingleCell"),
    .inflectionPoint
)

#' @rdname inflectionPoint
#' @export
setMethod(
    "inflectionPoint",
    signature("seurat"),
    .inflectionPoint
)
