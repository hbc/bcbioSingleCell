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
#' plotUMIsVsGenes(indrops_small)
NULL



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups = NULL,
        trendline = FALSE,
        color = getOption("bcbio.discrete.color", NULL),
        trans = "log2",
        title = "UMIs vs. genes"
    ) {
        # FIXME Use do.call
        .plotQCScatterplot(
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
    }
)
