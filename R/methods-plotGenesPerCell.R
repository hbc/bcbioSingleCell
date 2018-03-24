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
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' plotGenesPerCell(bcb_small)
#'
#' # seurat ====
#' plotGenesPerCell(seurat_small)
NULL



# Constructors =================================================================
.plotGenesPerCell <- function(
    object,
    geom = c("violin", "boxplot", "histogram", "ridgeline"),
    interestingGroups,
    min,
    max,
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



# Methods ======================================================================
#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("bcbioSingleCell"),
    .plotGenesPerCell
)



#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("seurat"),
    .plotGenesPerCell
)
