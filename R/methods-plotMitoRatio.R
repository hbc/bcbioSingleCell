#' Plot Mitochondrial Transcript Abundance
#'
#' @rdname plotMitoRatio
#' @name plotMitoRatio
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#'
#' @examples
#' # bcbioSingleCell
#' \dontrun{
#' plotMitoRatio(bcb)
#' }
#'
#' # data.frame
#' # This is recommended for advanced users only
#' \dontrun{
#' metrics <- metrics(bcb)
#' plotMitoRatio(metrics)
#' }
NULL



# Constructors ====
#' @importFrom viridis scale_fill_viridis
.plotMitoRatio <- function(
    object,
    geom = "violin",
    max = Inf,
    interestingGroups,
    multiplexed = FALSE,
    samplesOnYAxis = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)) {
    .checkGeom(geom)
    col <- "mitoRatio"
    min <- 0
    if (geom == "boxplot") {
        p <- .plotQCBoxplot(
            metrics = object,
            col = col,
            min = min,
            max  = max)
    } else if (geom == "histogram") {
        p <- .plotQCHistogram(
            metrics = object,
            col = col,
            min = min,
            max  = max)
    } else if (geom == "ridgeline") {
        p <- .plotQCRidgeline(
            metrics = object,
            col = col,
            min = min,
            max  = max)
    } else if (geom == "violin") {
        p <- .plotQCViolin(
            metrics = object,
            col = col,
            min = min,
            max = max)
    }

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

    # Facets
    facets <- NULL
    if (isTRUE(multiplexed) & length(unique(object[["description"]])) > 1) {
        facets <- c(facets, "description")
    }
    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(object)) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (!is.null(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    } else {
        # Add median labels
        if (geom != "histogram") {
            p <- p + .medianLabels(object, medianCol = col, digits = 2)
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
    function(
        object,
        geom = "violin",
        max,
        interestingGroups,
        filterCells = TRUE,
        aggregateReplicates = TRUE,
        samplesOnYAxis = TRUE,
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- metadata(object)[["interestingGroups"]]
        }
        if (missing(max)) {
            max <- metadata(object) %>%
                .[["filterParams"]] %>%
                .[["maxMitoRatio"]]
            if (is.null(max)) {
                max <- Inf
            }
        }
        multiplexed <- metadata(object)[["multiplexedFASTQ"]]
        metrics <- metrics(
            object,
            interestingGroups = interestingGroups,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates)
        .plotMitoRatio(
            object = metrics,
            geom = geom,
            max = max,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill,
            multiplexed = multiplexed)
    })



#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("data.frame"),
    .plotMitoRatio)



#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("seurat"),
    function(
        object,
        geom = "violin",
        max = Inf,
        interestingGroups,
        multiplexed = FALSE,
        samplesOnYAxis = TRUE,
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- slot(object, "misc") %>%
                .[["bcbio"]] %>%
                .[["interestingGroups"]]
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        .plotMitoRatio(
            object = metrics,
            geom = geom,
            max = max,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill,
            multiplexed = multiplexed)
    })
