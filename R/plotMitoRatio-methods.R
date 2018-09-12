#' Plot Mitochondrial Transcript Abundance
#'
#' @name plotMitoRatio
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#'
#' @return `ggplot`.
#'
#' @examples
#' plotMitoRatio(indrops_small)
NULL



#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("SingleCellExperiment"),
    function(
        object,
        geom = c("violin", "ridgeline", "ecdf", "histogram", "boxplot"),
        interestingGroups = NULL,
        max = 1L,
        fill = getOption("bcbio.discrete.fill", NULL),
        trans = "sqrt",
        title = "mito ratio"
    ) {
        geom <- match.arg(geom)
        assert_all_are_in_left_open_range(max, lower = 0L, upper = 1L)
        # FIXME Use do.call
        .plotQCMetric(
            object = object,
            metricCol = "mitoRatio",
            geom = geom,
            interestingGroups = interestingGroups,
            max = max,
            trans = trans,
            ratio = TRUE,
            fill = fill,
            title = title
        )
    }
)
