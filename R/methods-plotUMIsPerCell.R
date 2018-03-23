#' Plot UMIs per Cell
#'
#' Plot the universal molecular identifiers (UMIs) per cell.
#'
#' @rdname plotUMIsPerCell
#' @name plotUMIsPerCell
#' @family Quality Control Metrics
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
NULL



# Constructors =================================================================
.plotUMIsPerCell <- function(
    object,
    geom = "violin",
    min = 0L,
    interestingGroups,
    flip = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    if (missing(min)) {
        min <- metadata(object)[["filterParams"]][["minUMIs"]]
    }
    .plotQCGeom(
        metrics = metrics(object, interestingGroups = interestingGroups),
        metricCol = "nUMI",
        geom = geom,
        min = min,
        interestingGroups = interestingGroups,
        flip = flip,
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
