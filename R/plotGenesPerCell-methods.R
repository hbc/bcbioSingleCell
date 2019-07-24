#' @name plotGenesPerCell
#' @author Michael Steinbaugh, Rory Kirchner
#' @include globals.R
#' @inherit bioverbs::plotGenesPerCell
#'
#' @inheritParams acidplots::params
#' @inheritParams basejump::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(indrops)
#' plotGenesPerCell(indrops)
NULL



#' @rdname plotGenesPerCell
#' @name plotGenesPerCell
#' @importFrom bioverbs plotGenesPerCell
#' @usage plotGenesPerCell(object, ...)
#' @export
NULL



## Updated 2019-07-24.
`plotGenesPerCell,bcbioSingleCell` <-  # nolint
    function(
        object,
        geom,
        interestingGroups = NULL,
        min = 0L,
        max = Inf,
        trans = "log2",
        fill,
        title = "genes per cell"
    ) {
        geom <- match.arg(geom)
        do.call(
            what = .plotQCMetric,
            args = list(
                object = object,
                metricCol = "nGene",
                geom = geom,
                interestingGroups = interestingGroups,
                min = min,
                max = max,
                trans = trans,
                fill = fill,
                title = title
            )
        )
    }

formals(`plotGenesPerCell,bcbioSingleCell`)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(`plotGenesPerCell,bcbioSingleCell`)[["geom"]] <- geom



#' @rdname plotGenesPerCell
#' @export
setMethod(
    f = "plotGenesPerCell",
    signature = signature("bcbioSingleCell"),
    definition = `plotGenesPerCell,bcbioSingleCell`
)
