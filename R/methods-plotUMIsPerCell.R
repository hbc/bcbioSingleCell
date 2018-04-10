#' Plot UMIs per Cell
#'
#' Plot the universal molecular identifiers (UMIs) per cell.
#'
#' @name plotUMIsPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#' @param point Label either the "`knee`" or "`inflection`" point.
#' @param pointsPerSample Calculate the barcode ranks individually per sample.
#'   Generally recommended unless there's a reason to use a single cutoff point
#'   across all samples in the analysis.
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
        min,
        point = c("knee", "inflection"),
        pointsPerSample = TRUE,
        trans = "log10",
        color = scale_color_viridis(discrete = TRUE),
        fill = scale_fill_viridis(discrete = TRUE)
    ) {
        geom <- match.arg(geom)
        point <- match.arg(point)
        assert_is_a_bool(pointsPerSample)

        p <- .plotQCMetric(
            object = object,
            metricCol = "nUMI",
            geom = geom,
            interestingGroups = interestingGroups,
            min = min,
            trans = trans,
            fill = fill
        )

        if (isTRUE(pointsPerSample)) {
            FUN <- .labelBarcodeRanksPerSample
        } else {
            FUN <- .labelBarcodeRanks
        }

        p <- FUN(
            p = p,
            object = object,
            geom = geom,
            point = point
        )

        p
    }
)
