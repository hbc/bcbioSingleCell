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
#' plotGenesPerCell(bcb_small)
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
        geom = c("violin", "boxplot", "histogram", "ridgeline"),
        interestingGroups,
        min = 0L,
        max = Inf,
        fill = scale_fill_viridis(discrete = TRUE)
    ) {
        geom <- match.arg(geom)
        .plotQCMetric(
            object = object,
            metricCol = "nGene",
            geom = geom,
            interestingGroups = interestingGroups,
            min = min,
            max = max,
            trans = "identity",
            fill = fill
        )
    }
)
