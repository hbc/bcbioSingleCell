#' Plot Novelty Score
#'
#' "Novelty" refers to log10 genes detected per count.
#'
#' @name plotNovelty
#' @family Quality Control Functions
#' @author Michael Steinbaugh
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotNovelty(bcb_small)
#'
#' # SingleCellExperiment ====
#' plotNovelty(cellranger_small)
#'
#' # seurat ====
#' # `object@meta.data` must contain `log10GenesPerUMI`
#' \dontrun{
#' plotNovelty(Seurat::pbmc_small)
#' }
NULL



# Methods ======================================================================
#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    signature("SingleCellExperiment"),
    function(
        object,
        geom = c("violin", "boxplot", "histogram", "ridgeline"),
        interestingGroups,
        min = 0L,
        fill = scale_fill_hue()
    ) {
        assert_all_are_in_right_open_range(min, lower = 0L, upper = 1L)
        geom <- match.arg(geom)
        .plotQCMetric(
            object = object,
            metricCol = "log10GenesPerUMI",
            geom = geom,
            interestingGroups = interestingGroups,
            min = min,
            trans = "sqrt",
            fill = fill
        )
    }
)



#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    signature("seurat"),
    getMethod("plotNovelty", "SingleCellExperiment")
)
