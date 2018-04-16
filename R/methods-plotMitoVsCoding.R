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
#' plotMitoVsCoding(bcb_small)
#'
#' # SingleCellExperiment ====
#' plotMitoVsCoding(cellranger_small)
#'
#' # seurat ====
#' # `object@meta.data` must contain `nCoding`, `mitoRatio`
#' \dontrun{
#' plotMitoVsCoding(Seurat::pbmc_small)
#' }
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
        color = scale_color_viridis(discrete = TRUE)
    ) {
        .plotQCScatterplot(
            object = object,
            xCol = "nCoding",
            yCol = "nMito",
            xTrans = "log2",
            yTrans = "log2"
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
