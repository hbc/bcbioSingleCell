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
#' @param barcodeRanks Calculate barcode ranks per sample to label knee or
#'   inflection points. If `TRUE`, this will overwrite the `interestingGroups`
#'   setting and always plot per sample.
#' @param labelPoint Label either the "`knee`" or "`inflection`" point.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotUMIsPerCell(indrops_small)
#'
#' # Label the barcode ranks
#' plotUMIsPerCell(indrops_small, barcodeRanks = TRUE)
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
        barcodeRanks = TRUE,
        labelPoint = c("knee", "inflection"),
        trans = "log10",
        color = scale_color_hue(),
        fill = scale_fill_hue()
    ) {
        geom <- match.arg(geom)
        assert_is_a_bool(barcodeRanks)
        assert_is_a_bool(ranksPerSample)
        labelPoint <- match.arg(labelPoint)

        # Override interestingGroups if barcodeRanks is enabled
        if (isTRUE(barcodeRanks)) {
            interestingGroups <- "sampleName"
        }

        p <- .plotQCMetric(
            object = object,
            metricCol = "nUMI",
            geom = geom,
            interestingGroups = interestingGroups,
            min = min,
            trans = trans,
            color = color,
            fill = fill
        )

        # Calculate barcode ranks and label inflection or knee point
        if (isTRUE(barcodeRanks)) {
            p <- .labelBarcodeRanksPerSample(
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
