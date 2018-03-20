#' Plot Mitochondrial Transcript Abundance
#'
#' @name plotMitoRatio
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#'
#' @examples
#' # bcbioSingleCell ====
#' plotMitoRatio(bcb_small)
#'
#' # seurat ====
#' plotMitoRatio(seurat_small)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase interestingGroups
#' @importFrom ggplot2 coord_flip facet_wrap labs
.plotMitoRatio <- function(
    object,
    geom = "violin",
    max = Inf,
    interestingGroups,
    samplesOnYAxis = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    if (missing(max)) {
        max <- metadata(object)[["filterParams"]][["maxMitoRatio"]]
    }

    metrics <- metrics(object, interestingGroups)
    metricCol <- "mitoRatio"

    p <- .plotQCGeom(
        metrics = metrics,
        geom = geom,
        metricCol = metricCol,
        max = max
    )

    # Label interesting groups
    if (!missing(interestingGroups)) {
        p <- p + labs(fill = paste(interestingGroups, collapse = ":\n"))
    } else {
        p <- p + labs(color = NULL, fill = NULL)
    }

    # Color palette
    if (!is.null(fill)) {
        p <- p + fill
    }

    # Median labels
    if (geom %in% validMedianGeom) {
        p <- p + .medianLabels(metrics, medianCol = metricCol, digits = 2L)
    }

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(metrics))) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    if (isTRUE(samplesOnYAxis) && geom %in% validQCGeomFlip) {
        p <- p + coord_flip()
    }

    p
}



# Methods ======================================================================
#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("bcbioSingleCell"),
    .plotMitoRatio
)



#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("seurat"),
    .plotMitoRatio
)
