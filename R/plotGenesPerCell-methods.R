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
#' plotGenesPerCell(indrops_small)
NULL



#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("SingleCellExperiment"),
    function(
        object,
        geom = c("violin", "ridgeline", "ecdf", "histogram", "boxplot"),
        interestingGroups,
        min = 0L,
        max = Inf,
        trans = "log2",
        fill = getOption("bcbio.discrete.fill", NULL),
        title = "genes per cell"
    ) {
        geom <- match.arg(geom)
        .plotQCMetric(
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
    }
)
