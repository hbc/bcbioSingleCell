#' @name plotGenesPerCell
#' @author Michael Steinbaugh, Rory Kirchner
#' @include globals.R
#' @inherit bioverbs::plotGenesPerCell
#' @inheritParams basejump::params
#' @examples
#' data(indrops)
#' plotGenesPerCell(indrops)
NULL



#' @importFrom bioverbs plotGenesPerCell
#' @aliases NULL
#' @export
bioverbs::plotGenesPerCell



plotGenesPerCell.SingleCellExperiment <-  # nolint
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

formals(plotGenesPerCell.SingleCellExperiment)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotGenesPerCell.SingleCellExperiment)[["geom"]] <- geom



#' @rdname plotGenesPerCell
#' @export
setMethod(
    f = "plotGenesPerCell",
    signature = signature("SingleCellExperiment"),
    definition = plotGenesPerCell.SingleCellExperiment
)
