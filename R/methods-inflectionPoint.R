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
#' inflectionPoint(pbmc_small)
NULL



# Methods ======================================================================
#' @rdname inflectionPoint
#' @export
setMethod(
    "inflectionPoint",
    signature("SingleCellExperiment"),
    function(object) {
        validObject(object)
        metrics <- metrics(object)
        assert_is_subset("nUMI", colnames(metrics))
        counts <- sort(metrics[["nUMI"]])

        # Trigonometric calculations
        xdata <- seq_along(counts) / length(counts)
        cs <- cumsum(counts / sum(counts))
        m <- max(cs) / length(counts)
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
        inflection <- counts[which.max(dists)]
        names(inflection) <- "nUMI"
        inflection
    }
)



#' @rdname inflectionPoint
#' @export
setMethod(
    "inflectionPoint",
    signature("seurat"),
    getMethod("inflectionPoint", "SingleCellExperiment")
)
