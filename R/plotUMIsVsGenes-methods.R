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
#' # SingleCellExperiment ====
#' plotUMIsVsGenes(cellranger_small)
NULL



# Methods ======================================================================
#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups,
        trendline = FALSE,
        color = getOption("bcbio.discrete.color", NULL),
        trans = "log2",
        title = "UMIs vs. genes"
    ) {
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



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("seurat"),
    getMethod("plotUMIsVsGenes", "SingleCellExperiment")
)
