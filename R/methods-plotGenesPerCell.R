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
#' @param geom Plot type. Supported formats: `boxplot`, `histogram`,
#'   `ridgeline`, and `violin`.
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
#' \dontrun{
#' plotGenesPerCell(bcb)
#' }
#'
#' # seurat
#' \dontrun{
#' plotGenesPerCell(seurat)
#' }
#'
#' # metrics data.frame
#' \dontrun{
#' metrics <- metrics(bcb)
#' plotGenesPerCell(metrics)
#' }
NULL



# Constructors ====
#' @importFrom viridis scale_fill_viridis
.plotGenesPerCell <- function(
    object,
    geom = "boxplot",
    min = 0,
    max = Inf,
    interestingGroups,
    multiplexed = FALSE,
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
            p <- p + .medianLabels(object, medianCol = metricCol)
        }
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
        geom = "boxplot",
        min,
        max,
        interestingGroups,
        filterCells = TRUE,
        aggregateReplicates = TRUE,
        samplesOnYAxis = TRUE,
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- metadata(object)[["interestingGroups"]]
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
        multiplexed <- metadata(object)[["multiplexedFASTQ"]]
        .plotGenesPerCell(
            object = metrics,
            geom = geom,
            min = min,
            max = max,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill,
            multiplexed = multiplexed)
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
        geom = "boxplot",
        min,
        max,
        interestingGroups,
        multiplexed = FALSE,
        samplesOnYAxis = TRUE,
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- slot(object, "misc") %>%
                .[["bcbio"]] %>%
                .[["interestingGroups"]]
        }
        if (missing(min)) {
            min <- slot(object, "misc") %>%
                .[["bcbio"]] %>%
                .[["filterParams"]] %>%
                .[["minGenes"]]
        }
        if (missing(max)) {
            max <- slot(object, "misc") %>%
                .[["bcbio"]] %>%
                .[["filterParams"]] %>%
                .[["maxGenes"]]
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        .plotGenesPerCell(
            object = metrics,
            geom = geom,
            min = min,
            max = max,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill,
            multiplexed = multiplexed)
    })
