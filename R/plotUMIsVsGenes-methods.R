#' @name plotUMIsVsGenes
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit bioverbs::plotUMIsVsGenes
#' @inheritParams minimalism::params
#' @inheritParams basejump::params
#' @examples
#' data(indrops)
#' plotUMIsVsGenes(indrops)
NULL



#' @importFrom bioverbs plotUMIsVsGenes
#' @aliases NULL
#' @export
bioverbs::plotUMIsVsGenes



plotUMIsVsGenes.SingleCellExperiment <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        trendline = FALSE,
        color,
        trans = "log2",
        title = "UMIs vs. genes"
    ) {
        do.call(
            what = .plotQCScatterplot,
            args = list(
                object = object,
                interestingGroups = interestingGroups,
                trendline = trendline,
                xCol = "nUMI",
                yCol = "nGene",
                color = color,
                xTrans = trans,
                yTrans = trans,
                title = title
            )
        )
    }

formals(plotUMIsVsGenes.SingleCellExperiment)[["color"]] <-
    formalsList[["color.discrete"]]



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    f = "plotUMIsVsGenes",
    signature = signature("SingleCellExperiment"),
    definition = plotUMIsVsGenes.SingleCellExperiment
)
