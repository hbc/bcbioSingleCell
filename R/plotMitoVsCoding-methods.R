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



.plotMitoVsCoding.SCE <-  # nolint
    function(
        object,
        interestingGroups = NULL,
        trendline = FALSE,
        color = getOption("basejump.discrete.color", NULL),
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



#' @rdname plotMitoVsCoding
#' @export
setMethod(
    f = "plotMitoVsCoding",
    signature = signature("SingleCellExperiment"),
    definition = .plotMitoVsCoding.SCE
)
