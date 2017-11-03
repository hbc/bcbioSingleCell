#' Plot Mitochondrial Transcript Abundance
#'
#' @rdname plotMitoRatio
#' @name plotMitoRatio
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#'
#' @param color *Only applies to scatterplot*. Desired ggplot color scale.
#'   Defaults to [viridis::scale_color_viridis()]. Must supply discrete values.
#'   When set to `NULL`, the default ggplot2 color palette will be used. If
#'   manual color definitions are desired, we recommend using
#'   [ggplot2::scale_color_manual()].
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
    color = scale_color_viridis(discrete = TRUE),
    fill = scale_fill_viridis(discrete = TRUE)) {
    metricCol <- "mitoRatio"
    if (geom == "scatterplot") {
        p <- .plotQCScatterplot(
            object,
            xCol = "nCoding",
            yCol = "nMito")
    } else {
        p <- .dynamicQCPlot(
            object,
            metricCol = metricCol,
            min = 0,
            max = max,
            geom = geom)
    }

    # Label interesting groups
    if (!missing(interestingGroups)) {
        p <- p +
            labs(color = paste(interestingGroups, collapse = ":\n"),
                 fill = paste(interestingGroups, collapse = ":\n"))
    } else {
        p <- p + labs(color = NULL, fill = NULL)
    }

    # Color palette
    if (!is.null(color)) {
        p <- p + color
    }
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
        if (geom %in% validMedianGeom) {
            p <- p + .medianLabels(object, medianCol = metricCol, digits = 2)
        }
    }

    if (isTRUE(samplesOnYAxis) & geom %in% validQCGeomFlip) {
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
