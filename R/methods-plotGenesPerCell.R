#' Plot Genes per Cell
#'
#' @details A violin plot is a comact display of a continuous distribution. It
#'   is a blend of [geom_boxplot()] and [geom_density()]: a violin plot is a
#'   mirrored density plot displayed in the same way as a boxplot.
#'
#' @rdname plotGenesPerCell
#' @name plotGenesPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams AllGenerics
#' @inheritParams metrics
#' @inheritParams plotQC
#'
#' @param min Recommended minimum value cutoff.
#' @param max Recommended maximum value cutoff.
#' @param fill Desired ggplot fill scale. Defaults to
#'   [viridis::scale_fill_viridis()]. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using [ggplot2::scale_fill_manual()].
#' @param samplesOnYAxis Plot the samples on the y axis. Doesn't apply to
#'   histogram.
#'
#' @return [ggplot].
#'
#' @examples
#' # bcbioSingleCell
#' bcb <- examples[["bcb"]]
#' plotGenesPerCell(bcb)
#'
#' # seurat
#' seurat <- examples[["seurat"]]
#' plotGenesPerCell(seurat)
#'
#' # data.frame
#' metrics <- metrics(bcb)
#' plotGenesPerCell(metrics)
NULL



# Constructors ====
#' @importFrom viridis scale_fill_viridis
.plotGenesPerCell <- function(
    object,
    geom = "violin",
    min = 0,
    max = Inf,
    interestingGroups,
    samplesOnYAxis = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)) {
    metricCol <- "nGene"
    p <- .plotQCGeom(
        object,
        geom = geom,
        metricCol = metricCol,
        min = min,
        max = max)

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
        p <- p + .medianLabels(object, medianCol = metricCol)
    }

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
        facets <- "sampleNameAggregate"
    }
    if (!is.null(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    if (isTRUE(samplesOnYAxis) & geom %in% validQCGeomFlip) {
        p <- p + coord_flip()
    }

    p
}



# Methods ====
#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("bcbioSingleCell"),
    function(
        object,
        geom = "violin",
        min,
        max,
        interestingGroups,
        samplesOnYAxis = TRUE,
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        if (missing(min)) {
            min <- metadata(object)[["filterParams"]][["minGenes"]]
        }
        if (missing(max)) {
            max <- metadata(object)[["filterParams"]][["maxGenes"]]
        }
        metrics <- metrics(
            object,
            interestingGroups = interestingGroups)
        .plotGenesPerCell(
            object = metrics,
            geom = geom,
            min = min,
            max = max,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill)
    })



#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("data.frame"),
    .plotGenesPerCell)



#' @rdname plotGenesPerCell
#' @export
setMethod(
    "plotGenesPerCell",
    signature("seurat"),
    function(
        object,
        geom = "violin",
        min,
        max,
        interestingGroups,
        samplesOnYAxis = TRUE,
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        if (missing(min)) {
            min <- bcbio(object)[["filterParams"]][["minGenes"]]
        }
        if (missing(max)) {
            max <- bcbio(object)[["filterParams"]][["maxGenes"]]
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        .plotGenesPerCell(
            object = metrics,
            geom = geom,
            min = min,
            max = max,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill)
    })
