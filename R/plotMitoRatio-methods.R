#' @name plotMitoRatio
#' @author Michael Steinbaugh, Rory Kirchner
#' @include globals.R
#' @inherit basejump::plotMitoRatio
#' @inheritParams basejump::params
#' @examples
#' data(indrops)
#' plotMitoRatio(indrops)
NULL



#' @importFrom basejump plotMitoRatio
#' @aliases NULL
#' @export
basejump::plotMitoRatio



plotMitoRatio.SingleCellExperiment <-  # nolint
    function(
        object,
        geom,
        interestingGroups = NULL,
        max = 1L,
        fill,
        trans = "sqrt",
        title = "mito ratio"
    ) {
        assert(isInLeftOpenRange(max, lower = 0L, upper = 1L))
        geom <- match.arg(geom)
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

formals(plotMitoRatio.SingleCellExperiment)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotMitoRatio.SingleCellExperiment)[["geom"]] <- geom



#' @rdname plotMitoRatio
#' @export
setMethod(
    f = "plotMitoRatio",
    signature = signature("SingleCellExperiment"),
    definition = plotMitoRatio.SingleCellExperiment
)
