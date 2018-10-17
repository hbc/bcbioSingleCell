#' Plot Genes per Cell
#'
#' @name plotGenesPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' data(indrops_small)
#' plotGenesPerCell(indrops_small)
NULL



.plotGenesPerCell.SCE <-  # nolint
    function(
        object,
        geom = c("violin", "ridgeline", "ecdf", "histogram", "boxplot"),
        interestingGroups = NULL,
        min = 0L,
        max = Inf,
        trans = "log2",
        fill = getOption("basejump.discrete.fill", NULL),
        title = "genes per cell"
    ) {
        geom <- match.arg(geom)
        do.call(
            what = .plotQCMetric,
            args = list(
                object = object,
                metricCol = "nGene",
                geom = geom,
                interestingGroups = interestingGroups,
                min = min,
                max = max,
                trans = trans,
                fill = fill,
                title = title
            )
        )
    }



#' @rdname plotGenesPerCell
#' @export
setMethod(
    f = "plotGenesPerCell",
    signature = signature("SingleCellExperiment"),
    definition = .plotGenesPerCell.SCE
)
