#' Plot UMI and Gene Correlation
#'
#' @name plotUMIsVsGenes
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' data(indrops_small)
#' plotUMIsVsGenes(indrops_small)
NULL



.plotUMIsVsGenes.SCE <-  # nolint
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



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    f = "plotUMIsVsGenes",
    signature = signature("SingleCellExperiment"),
    definition = .plotUMIsVsGenes.SCE
)
