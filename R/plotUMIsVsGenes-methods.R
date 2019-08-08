#' @name plotUMIsVsGenes
#' @author Michael Steinbaugh, Rory Kirchner
#' @inherit bioverbs::plotUMIsVsGenes
#' @note Updated 2019-08-08.
#'
#' @inheritParams acidroxygen::params
#' @param ... Additional arguments.
#'
#' @examples
#' data(indrops)
#' plotUMIsVsGenes(indrops)
NULL



#' @rdname plotUMIsVsGenes
#' @name plotUMIsVsGenes
#' @importFrom bioverbs plotUMIsVsGenes
#' @usage plotUMIsVsGenes(object, ...)
#' @export
NULL



## Updated 2019-07-24.
`plotUMIsVsGenes,bcbioSingleCell` <-  # nolint
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
                ## Note that these were renamed in v0.3.19 to better match
                ## conventions used in Chromium and Seurat packages.
                xCol = "nCount",
                yCol = "nFeature",
                color = color,
                xTrans = trans,
                yTrans = trans,
                title = title
            )
        )
    }

formals(`plotUMIsVsGenes,bcbioSingleCell`)[["color"]] <-
    formalsList[["color.discrete"]]



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    f = "plotUMIsVsGenes",
    signature = signature("bcbioSingleCell"),
    definition = `plotUMIsVsGenes,bcbioSingleCell`
)
