#' Calculate the Cellular Barcode Read Depth Inflection Point
#'
#' @name inflectionPoint
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
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
#' inflectionPoint(Seurat::pbmc_small)
NULL



# Methods ======================================================================
#' @rdname inflectionPoint
#' @export
setMethod(
    "inflectionPoint",
    signature("SingleCellExperiment"),
    function(object) {
        .Deprecated("barcodeRanks")

        validObject(object)
        counts <- counts(object)

        # Calculate the total number of UMIs per cell (nUMI)
        totals <- colSums(counts) %>%
            sort(decreasing = FALSE)

        # Trigonometric calculations
        xdata <- seq_along(totals) / length(totals)
        cs <- cumsum(totals / sum(totals))
        m <- max(cs) / length(totals)
        dists <- mcmapply(
            xdata = xdata,
            cs = cs,
            MoreArgs = list(m = m),
            FUN = function(xdata, cs, m) {
                sin(atan(m)) * abs(cs - xdata)
            },
            SIMPLIFY = TRUE,
            USE.NAMES = FALSE
        )

        # Return the inflection point as the expression value
        totals[which.max(dists)]
    }
)



#' @rdname inflectionPoint
#' @export
setMethod(
    "inflectionPoint",
    signature("seurat"),
    getMethod("inflectionPoint", "SingleCellExperiment")
)
