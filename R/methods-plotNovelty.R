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
        geom = c("violin", "boxplot", "histogram", "ridgeline"),
        interestingGroups,
        min,
        fill = scale_fill_viridis(discrete = TRUE)
    ) {
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
