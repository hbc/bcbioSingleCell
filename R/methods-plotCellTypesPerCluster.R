#' Plot Cell Types per Cluster
#'
#' Plot the geometric mean of the significant marker genes for every known cell
#' type (per unbiased cluster). Cell types with too few (`min` cutoff) or too
#' many (`max` cutoff) marker genes will be skipped.
#'
#' @rdname plotCellTypesPerCluster
#' @name plotCellTypesPerCluster
#' @author Michael Steinbaugh
#'
#' @inheritParams AllGenerics
#'
#' @param min Minimum number of marker genes per cell type.
#' @param max Maximum number of marker genes per cell type.
#'
#' @return [ggplot].
#'
#' @examples
#' \dontrun{
#' cells <- cellTypesPerCluster(knownMakersDetected)
#' plotCellTypesPerCluster(cells)
#' }
NULL



# Constructors ====
.plotCellTypesPerCluster <- function(object) {
    object
    plotMarkerTSNE()
}



# Methods ====
#' @rdname plotCellTypesPerCluster
#' @export
