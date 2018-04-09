#' Plot Mitochondrial Transcript Abundance
#'
#' @name plotMitoRatio
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotMitoRatio(bcb_small)
#'
#' # SingleCellExperiment ====
#' plotMitoRatio(cellranger_small)
#'
#' # seurat ====
#' plotMitoRatio(seurat_small)
NULL



# Methods ======================================================================
#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("SingleCellExperiment"),
    function(
        object,
        geom = c("violin", "boxplot", "histogram", "ridgeline"),
        interestingGroups,
        max,
        fill = scale_fill_viridis(discrete = TRUE)
    ) {
        geom <- match.arg(geom)
        .plotQCMetric(
            object = object,
            metricCol = "mitoRatio",
            geom = geom,
            interestingGroups = interestingGroups,
            max = max,
            trans = "identity",
            fill = fill
        )
    }
)
