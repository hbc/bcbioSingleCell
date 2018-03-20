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
#' # bcbioSingleCell ====
#' plotUMIsPerCell(bcb_small)
#'
#' # seurat ====
#' plotUMIsPerCell(seurat_small)
NULL



# Constructors =================================================================
#' @importFrom ggplot2 coord_flip facet_wrap
.plotUMIsPerCell <- function(
    object,
    geom = "violin",
    min = 0L,
    interestingGroups,
    samplesOnYAxis = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    if (missing(min)) {
        min <- metadata(object)[["filterParams"]][["minUMIs"]]
    }

    metrics <- metrics(object, interestingGroups)
    metricCol <- "nUMI"

    p <- .plotQCGeom(
        metrics = metrics,
        geom = geom,
        metricCol = metricCol,
        min = min
    )

    # Label interesting groups
    if (!missing(interestingGroups)) {
        p <- p + labs(fill = paste(interestingGroups, collapse = ":\n"))
    } else {
        p <- p + labs(fill = NULL)
    }

    # Color palette
    if (!is.null(fill)) {
        p <- p + fill
    }

    # Median labels
    if (geom %in% validMedianGeom) {
        p <- p + .medianLabels(metrics, medianCol = metricCol)
    }

    # Wrap aggregated samples
    facets <- NULL
    if (isTRUE(.checkAggregate(metrics))) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    # Flip axis, if desired
    if (isTRUE(samplesOnYAxis) && geom %in% validQCGeomFlip) {
        p <- p + coord_flip()
    }

    p
}



# Methods ======================================================================
#' @rdname plotUMIsPerCell
#' @importFrom viridis scale_fill_viridis
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
