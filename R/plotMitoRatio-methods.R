#' @name plotMitoRatio
#' @author Michael Steinbaugh, Rory Kirchner
#' @include globals.R
#' @inherit bioverbs::plotMitoRatio
#'
#' @inheritParams acidplots::params
#' @inheritParams basejump::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(indrops)
#' plotMitoRatio(indrops)
NULL



#' @rdname plotMitoRatio
#' @name plotMitoRatio
#' @importFrom bioverbs plotMitoRatio
#' @usage plotMitoRatio(object, ...)
#' @export
NULL



## Updated 2019-07-24.
`plotMitoRatio,bcbioSingleCell` <-  # nolint
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

formals(`plotMitoRatio,bcbioSingleCell`)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(`plotMitoRatio,bcbioSingleCell`)[["geom"]] <- geom



#' @rdname plotMitoRatio
#' @export
setMethod(
    f = "plotMitoRatio",
    signature = signature("bcbioSingleCell"),
    definition = `plotMitoRatio,bcbioSingleCell`
)
