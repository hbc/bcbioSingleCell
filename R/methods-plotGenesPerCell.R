#' Plot Genes per Cell
#'
#' @rdname plotGenesPerCell
#' @name plotGenesPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inheritParams AllGenerics
#' @inheritParams metrics
#'
#' @param geom Plot type. Currently support formats: `boxplot`, `histogram`,
#'   `ridgeline`, and `violin` (default).
#' @param min Recommended minimum value cutoff.
#' @param max Recommended maximum value cutoff.
#' @param fill Desired ggplot fill scale. Defaults to
#'   [viridis::scale_fill_viridis()]. Must supply discrete values. When set to
#'   `NULL`, the default ggplot2 color palette will be used. If manual color
#'   definitions are desired, we recommend using [ggplot2::scale_fill_manual()].
#' @param samplesOnYAxis Plot the samples on the y axis. Doesn't apply to
#'   histogram.
#'
#' @return [ggplot] grid.
NULL



# Constructors ====
#' @importFrom viridis scale_fill_viridis
.plotGenesPerCell <- function(
    object,
    interestingGroups,
    geom = "violin",
    min,
    max,
    filterCells = TRUE,
    aggregateReplicates = TRUE,
    samplesOnYAxis = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)) {
    .checkGeom(geom)
    if (missing(interestingGroups)) {
        interestingGroups <-
            metadata(object)[["interestingGroups"]]
    }
    if (missing(min)) {
        min <- metadata(object) %>%
            .[["filterParams"]] %>%
            .[["minGenes"]]
        if (is.null(min)) {
            min <- 0
        }
    }
    if (missing(max)) {
        max <- metadata(object) %>%
            .[["filterParams"]] %>%
            .[["maxGenes"]]
        if (is.null(max)) {
            max <- Inf
        }
    }
    metrics <- metrics(
        object,
        interestingGroups = interestingGroups,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates)
    col <- "nGene"
    if (geom == "boxplot") {
        p <- .plotQCBoxplot(
            metrics = metrics,
            col = col,
            min = min,
            max  = max)
    } else if (geom == "histogram") {
        p <- .plotQCHistogram(
            metrics = metrics,
            col = col,
            min = min,
            max  = max)
    } else if (geom == "ridgeline") {
        p <- .plotQCRidgeline(
            metrics = metrics,
            col = col,
            min = min,
            max  = max)
    } else if (geom == "violin") {
        p <- .plotQCViolin(
            metrics = metrics,
            col = col,
            min = min,
            max = max)
    }

    # Label interesting groups
    p <- p + labs(fill = paste(interestingGroups, collapse = ":\n"))

    # Color palette
    if (!is.null(fill)) {
        p <- p + fill
    }

    # Facets
    facets <- NULL
    if (isTRUE(metadata(object)[["multiplexedFASTQ"]]) &
        length(unique(metrics[["description"]])) > 1) {
        facets <- c(facets, "description")
    }
    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(metrics)) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (!is.null(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    } else {
        # Add median labels
        if (geom != "histogram") {
            p <- p + .medianLabels(metrics, medianCol = col, digits = 0)
        }
    }

    if (isTRUE(samplesOnYAxis) & geom %in% validGCGeomFlip) {
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
    .plotGenesPerCell)
