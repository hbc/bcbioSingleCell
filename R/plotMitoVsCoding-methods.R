#' Plot Mitochondrial vs. Coding Counts
#'
#' @name plotMitoVsCoding
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotMitoVsCoding(indrops_small)
NULL



# Methods ======================================================================
#' @rdname plotMitoVsCoding
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("SingleCellExperiment"),
    function(
        object,
        interestingGroups,
        trendline = FALSE,
        color = getOption("bcbio.discrete.color", NULL),
        trans = "log2",
        title = "mito vs. coding"
    ) {
        .plotQCScatterplot(
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
    }
)
