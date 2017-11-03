#' Plot Mitochondrial Transcript Abundance
#'
#' @rdname plotMitoRatio
#' @name plotMitoRatio
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
NULL



# Constructors ====
#' @importFrom viridis scale_fill_viridis
.plotMitoRatio <- function(
    object,
    interestingGroups,
    geom = "violin",
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
    if (missing(max)) {
        max <- metadata(object) %>%
            .[["filterParams"]] %>%
            .[["maxMitoRatio"]]
        if (is.null(max)) {
            max <- Inf
        }
    }
    min <- 0
    metrics <- metrics(
        object,
        filterCells = filterCells,
        aggregateReplicates = aggregateReplicates) %>%
        uniteInterestingGroups(interestingGroups)
    col <- "mitoRatio"
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
            p <- p + .medianLabels(metrics, medianCol = col, digits = 2)
        }
    }

    if (isTRUE(samplesOnYAxis) & geom %in% validGCGeomFlip) {
        p <- p + coord_flip()
    }

    p
}



# Methods ====
#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("bcbioSingleCell"),
    .plotMitoRatio)
