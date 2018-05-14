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
#' # bcbioSingleCell ====
#' plotGenesPerCell(indrops_small)
#'
#' # SingleCellExperiment ====
#' plotGenesPerCell(cellranger_small)
#'
#' # seurat ====
#' plotGenesPerCell(Seurat::pbmc_small)
NULL



# Methods ======================================================================
#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("SingleCellExperiment"),
    function(
        object,
        geom = c("ecdf", "ridgeline", "violin", "histogram", "boxplot"),
        interestingGroups,
        min = 0L,
        max = Inf,
        trans = "log2",
        fill = scale_fill_hue(),
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



#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("seurat"),
    getMethod("plotGenesPerCell", "SingleCellExperiment")
)
