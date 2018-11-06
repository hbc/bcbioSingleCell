#' Plot Mitochondrial Transcript Abundance
#'
#' @name plotMitoRatio
#' @author Michael Steinbaugh, Rory Kirchner
#' @include globals.R
#'
#' @inheritParams basejump::params
#'
#' @return `ggplot`.
#'
#' @examples
#' data(indrops_small)
#' plotMitoRatio(indrops_small)
NULL



plotMitoRatio.SingleCellExperiment <-  # nolint
    function(
        object,
        geom,
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
formals(plotMitoRatio.SingleCellExperiment)[["geom"]] <- geom



#' @rdname plotMitoRatio
#' @export
setMethod(
    f = "plotMitoRatio",
    signature = signature("SingleCellExperiment"),
    definition = plotMitoRatio.SingleCellExperiment
)
