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
#' # `object@meta.data` must contain `mitoRatio`
#' \dontrun{
#' plotMitoRatio(Seurat::pbmc_small)
#' }
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
        max = 1L,
        fill = scale_fill_viridis(discrete = TRUE)
    ) {
        geom <- match.arg(geom)
        assert_all_are_in_left_open_range(max, lower = 0L, upper = 1L)
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



#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("seurat"),
    getMethod("plotMitoRatio", "SingleCellExperiment")
)
