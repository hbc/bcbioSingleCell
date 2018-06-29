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
#' # SingleCellExperiment ====
#' plotMitoRatio(cellranger_small)
NULL



# Methods ======================================================================
#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("SingleCellExperiment"),
    function(
        object,
        geom = c("ecdf", "ridgeline", "violin", "histogram", "boxplot"),
        interestingGroups,
        max = 1L,
        fill = NULL,
        trans = "sqrt",
        title = "mito ratio"
    ) {
        geom <- match.arg(geom)
        assert_all_are_in_left_open_range(max, lower = 0L, upper = 1L)
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



#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("seurat"),
    getMethod("plotMitoRatio", "SingleCellExperiment")
)
