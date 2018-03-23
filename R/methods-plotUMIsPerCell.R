#' Plot UMIs per Cell
#'
#' Plot the universal molecular identifiers (UMIs) per cell.
#'
#' @name plotUMIsPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#'
#' @examples
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' plotUMIsPerCell(bcb_small)
#'
#' # seurat ====
#' plotUMIsPerCell(seurat_small)
NULL



# Constructors =================================================================
.plotUMIsPerCell <- function(
    object,
    geom = c("violin", "boxplot", "histogram", "ridgeline"),
    interestingGroups,
    min,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    geom <- match.arg(geom)
    .plotQCMetric(
        object = object,
        metricCol = "nUMI",
        geom = geom,
        interestingGroups = interestingGroups,
        min = min,
        fill = fill
    )
}



# Methods ======================================================================
#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("bcbioSingleCell"),
    .plotUMIsPerCell
)



#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("seurat"),
    .plotUMIsPerCell
)
