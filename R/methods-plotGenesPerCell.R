#' Plot Genes per Cell
#'
#' @details A violin plot is a comact display of a continuous distribution. It
#'   is a blend of [geom_boxplot()] and [geom_density()]: a violin plot is a
#'   mirrored density plot displayed in the same way as a boxplot.
#'
#' @name plotGenesPerCell
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams general
#' @inheritParams metrics
#' @inheritParams plotQC
#' @param min Recommended minimum value cutoff.
#' @param max Recommended maximum value cutoff.
#' @param fill Desired ggplot fill scale. Defaults to
#'   [viridis::scale_fill_viridis()]. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using [ggplot2::scale_fill_manual()].
#' @param samplesOnYAxis Plot the samples on the y axis. Doesn't apply to
#'   histogram.
#'
#' @return `ggplot`.
#'
#' @examples
#' # bcbioSingleCell ====
#' plotGenesPerCell(bcb_small)
#'
#' # seurat ====
#' plotGenesPerCell(seurat_small)
NULL



# Constructors =================================================================
.plotGenesPerCell <- function(
    object,
    geom = "violin",
    min = 0L,
    max = Inf,
    interestingGroups,
    samplesOnYAxis = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    assert_is_a_bool(samplesOnYAxis)
    assertIsColorScaleDiscreteOrNULL(fill)

    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    if (missing(min)) {
        min <- metadata(object)[["filterParams"]][["minGenes"]]
    }
    if (missing(max)) {
        max <- metadata(object)[["filterParams"]][["maxGenes"]]
    }

    metrics <- metrics(object, interestingGroups)
    metricCol <- "nGene"

    p <- .plotQCGeom(
        metrics = metrics,
        geom = geom,
        metricCol = metricCol,
        min = min,
        max = max
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

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(metrics))) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    if (isTRUE(samplesOnYAxis) & geom %in% validQCGeomFlip) {
        p <- p + coord_flip()
    }

    p
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
