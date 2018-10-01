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



.plotNovelty.SCE <-  # nolint
    function(
        object,
        geom = c("violin", "ridgeline", "ecdf", "histogram", "boxplot"),
        interestingGroups = NULL,
        min = 0L,
        fill = getOption("basejump.discrete.fill", NULL),
        trans = "identity",
        title = "novelty : genes per UMI"
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



#' @rdname plotNovelty
#' @export
setMethod(
    f = "plotNovelty",
    signature = signature("SingleCellExperiment"),
    definition = .plotNovelty.SCE
)
