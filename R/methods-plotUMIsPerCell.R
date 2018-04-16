#' Plot UMIs per Cell
#'
#' Plot the universal molecular identifiers (UMIs) per cell.
#'
#' @name plotUMIsPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#'
#' @inheritParams general
#' @param barcodeRanks Calculate barcode ranks to label knee or inflection
#'   points.
#' @param ranksPerSample Calculate the barcode ranks individually per sample.
#'   Generally recommended unless there's a reason to use a single cutoff point
#'   across all samples in the analysis.
#' @param labelPoint Label either the "`knee`" or "`inflection`" point.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotUMIsPerCell(bcb_small)
#'
#' # SingleCellExperiment ====
#' plotUMIsPerCell(cellranger_small)
#'
#' # seurat ====
#' plotUMIsPerCell(Seurat::pbmc_small)
NULL



# Methods ======================================================================
#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("SingleCellExperiment"),
    function(
        object,
        geom = c("ecdf", "histogram", "ridgeline", "boxplot", "violin"),
        interestingGroups,
        min = 0L,
        barcodeRanks = FALSE,
        ranksPerSample = TRUE,
        labelPoint = c("knee", "inflection"),
        trans = "log10",
        color = scale_color_viridis(discrete = TRUE),
        fill = scale_fill_viridis(discrete = TRUE)
    ) {
        geom <- match.arg(geom)
        assert_is_a_bool(barcodeRanks)
        assert_is_a_bool(ranksPerSample)
        labelPoint <- match.arg(labelPoint)

        p <- .plotQCMetric(
            object = object,
            metricCol = "nUMI",
            geom = geom,
            interestingGroups = interestingGroups,
            min = min,
            trans = trans,
            fill = fill
        )

        # Calculate barcode ranks and label inflection or knee point
        if (isTRUE(barcodeRanks)) {
            if (isTRUE(ranksPerSample)) {
                fun <- .labelBarcodeRanksPerSample
            } else {
                fun <- .labelBarcodeRanks
            }
            p <- fun(
                p = p,
                object = object,
                geom = geom,
                point = labelPoint
            )
        }

        p
    }
)



#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("seurat"),
    getMethod("plotUMIsPerCell", "SingleCellExperiment")
)
