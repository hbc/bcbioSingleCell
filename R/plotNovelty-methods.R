#' @name plotNovelty
#' @author Michael Steinbaugh
#' @include globals.R
#' @inherit bioverbs::plotNovelty
#' @inheritParams basejump::params
#' @examples
#' data(indrops)
#' plotNovelty(indrops)
NULL



#' @importFrom bioverbs plotNovelty
#' @aliases NULL
#' @export
bioverbs::plotNovelty



plotNovelty.SingleCellExperiment <-  # nolint
    function(
        object,
        geom,
        interestingGroups = NULL,
        min = 0L,
        fill,
        trans = "identity",
        title = "novelty : genes per UMI"
    ) {
        assert(isInRightOpenRange(min, lower = 0L, upper = 1L))
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

formals(plotNovelty.SingleCellExperiment)[["fill"]] <-
    formalsList[["fill.discrete"]]
formals(plotNovelty.SingleCellExperiment)[["geom"]] <- geom



#' @rdname plotNovelty
#' @export
setMethod(
    f = "plotNovelty",
    signature = signature("SingleCellExperiment"),
    definition = plotNovelty.SingleCellExperiment
)
