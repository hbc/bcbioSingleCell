#' @name plotMitoVsCoding
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit bioverbs::plotMitoVsCoding
#' @inheritParams basejump::params
#' @examples
#' data(indrops)
#' plotMitoVsCoding(indrops)
NULL



#' @importFrom bioverbs plotMitoVsCoding
#' @aliases NULL
#' @export
bioverbs::plotMitoVsCoding



plotMitoVsCoding.SingleCellExperiment <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        trendline = FALSE,
        color,
        trans = "log2",
        title = "mito vs. coding"
    ) {
        do.call(
            what = .plotQCScatterplot,
            args = list(
                object = object,
                interestingGroups = interestingGroups,
                trendline = trendline,
                xCol = "nCoding",
                yCol = "nMito",
                color = color,
                xTrans = trans,
                yTrans = trans,
                title = title
            )
        )
    }

formals(plotMitoVsCoding.SingleCellExperiment)[["color"]] <-
    formalsList[["color.discrete"]]



#' @rdname plotMitoVsCoding
#' @export
setMethod(
    f = "plotMitoVsCoding",
    signature = signature("SingleCellExperiment"),
    definition = plotMitoVsCoding.SingleCellExperiment
)
