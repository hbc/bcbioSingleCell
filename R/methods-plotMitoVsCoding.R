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
#' # bcbioSingleCell ====
#' plotMitoVsCoding(indrops_small)
#'
#' # SingleCellExperiment ====
#' plotMitoVsCoding(cellranger_small)
#'
#' # seurat ====
#' plotMitoVsCoding(seurat_small)
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
        color = scale_color_hue(),
        trans = "log2",
        title = "mito vs. coding"
    ) {
        .plotQCScatterplot(
            object = object,
            interestingGroups = interestingGroups,
            trendline = trendline,
            xCol = "nCoding",
            yCol = "nMito",
            xTrans = trans,
            yTrans = trans,
            title = title
        )
    }
)



#' @rdname plotMitoVsCoding
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("seurat"),
    getMethod("plotMitoVsCoding", "SingleCellExperiment")
)
