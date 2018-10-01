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
#' plotMitoRatio(indrops_small)
NULL



.plotMitoRatio.SCE <-  # nolint
    function(
        object,
        geom = c("violin", "ridgeline", "ecdf", "histogram", "boxplot"),
        interestingGroups = NULL,
        max = 1L,
        fill = getOption("basejump.discrete.fill", NULL),
        trans = "sqrt",
        title = "mito ratio"
    ) {
        geom <- match.arg(geom)
        assert_all_are_in_left_open_range(max, lower = 0L, upper = 1L)
        do.call(
            what = .plotQCMetric,
            args = list(
                object = object,
                metricCol = "mitoRatio",
                geom = geom,
                interestingGroups = interestingGroups,
                max = max,
                trans = trans,
                ratio = TRUE,
                fill = fill,
                title = title
            )
        )
    }



#' @rdname plotMitoRatio
#' @export
setMethod(
    f = "plotMitoRatio",
    signature = signature("SingleCellExperiment"),
    definition = .plotMitoRatio.SCE
)
