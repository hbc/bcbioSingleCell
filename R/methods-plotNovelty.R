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
#' plotNovelty(indrops_small)
#'
#' # SingleCellExperiment ====
#' plotNovelty(cellranger_small)
#'
#' # seurat ====
#' plotNovelty(seurat_small)
NULL



# Methods ======================================================================
#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    signature("SingleCellExperiment"),
    function(
        object,
        geom = c("ecdf", "histogram", "ridgeline", "violin", "boxplot"),
        interestingGroups,
        min = 0L,
        fill = scale_fill_hue(),
        trans = "sqrt",
        title = "genes per UMI (novelty)"
    ) {
        assert_all_are_in_right_open_range(min, lower = 0L, upper = 1L)
        geom <- match.arg(geom)
        .plotQCMetric(
            object = object,
            metricCol = "log10GenesPerUMI",
            geom = geom,
            interestingGroups = interestingGroups,
            min = min,
            max = 1L,
            trans = trans,
            ratio = TRUE,
            fill = fill,
            title = title
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
