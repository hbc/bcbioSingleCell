#' Plot UMIs per Cell
#'
#' Plot the universal molecular identifiers (UMIs) per cell.
#'
#' @name plotUMIsPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#'
#' @examples
#' # bcbioSingleCell ====
#' plotUMIsPerCell(bcb_small)
#'
#' # seurat ====
#' plotUMIsPerCell(seurat_small)
#' plotUMIsPerCell(Seurat::pbmc_small)
NULL



# Constructors =================================================================
.plotUMIsPerCell <- function(
    object,
    geom = c("ecdf", "histogram", "boxplot", "ridgeline", "violin"),
    interestingGroups,
    min,
    trans = "log10",
    color = scale_color_viridis(discrete = TRUE),
    fill = scale_fill_viridis(discrete = TRUE)
) {
    geom <- match.arg(geom)
    .plotQCMetric(
        object = object,
        metricCol = "nUMI",
        geom = geom,
        interestingGroups = interestingGroups,
        min = min,
        trans = trans,
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
