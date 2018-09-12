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
#' plotNovelty(indrops_small)
NULL



#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    signature("SingleCellExperiment"),
    function(
        object,
        geom = c("violin", "ridgeline", "ecdf", "histogram", "boxplot"),
        interestingGroups = NULL,
        min = 0L,
        fill = getOption("bcbio.discrete.fill", NULL),
        trans = "identity",
        title = "genes per UMI (novelty)"
    ) {
        assert_all_are_in_right_open_range(min, lower = 0L, upper = 1L)
        geom <- match.arg(geom)
        do.call(
            what = .plotQCMetric,
            args = list(
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
        )
    }
)
