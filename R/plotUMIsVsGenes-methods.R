#' @name plotUMIsVsGenes
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit basejump::plotUMIsVsGenes
#' @inheritParams basejump::params
#' @examples
#' data(indrops)
#' plotUMIsVsGenes(indrops)
NULL



#' @importFrom basejump plotUMIsVsGenes
#' @aliases NULL
#' @export
basejump::plotUMIsVsGenes



# TODO Define the color globally.
plotUMIsVsGenes.SingleCellExperiment <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        trendline = FALSE,
        color = getOption("basejump.discrete.color", NULL),
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

# FIXME Need to set the formals here.



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    f = "plotUMIsVsGenes",
    signature = signature("SingleCellExperiment"),
    definition = plotUMIsVsGenes.SingleCellExperiment
)
