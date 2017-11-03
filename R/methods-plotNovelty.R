#' Plot Novelty Score
#'
#' "Novelty" refers to log10 genes detected per count.
#'
#' @rdname plotNovelty
#' @name plotNovelty
#' @family Quality Control Metrics
#' @author Michael Steinbaugh
#'
#' @inherit plotGenesPerCell
#'
#' @examples
#' # bcbioSingleCell
#' \dontrun{
#' plotNovelty(bcb)
#' }
#'
#' # seurat
#' \dontrun{
#' plotNovelty(seurat)
#' }
#'
#' # data.frame
#' \dontrun{
#' metrics <- metrics(bcb)
#' plotNovelty(metrics)
#' }
NULL



# Constructors ====
#' @importFrom viridis scale_fill_viridis
.plotNovelty <- function(
    object,
    geom = "boxplot",
    min = 0,
    interestingGroups,
    multiplexed = FALSE,
    samplesOnYAxis = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)) {
    metricCol <- "log10GenesPerUMI"
    p <- .plotQCGeom(
        object,
        geom = geom,
        metricCol = metricCol,
        min = min)

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
#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    signature("bcbioSingleCell"),
    function(
        object,
        geom = "boxplot",
        min,
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
                .[["minNovelty"]]
            if (is.null(min)) {
                min <- 0
            }
        }
        multiplexed <- metadata(object)[["multiplexedFASTQ"]]
        metrics <- metrics(
            object,
            interestingGroups = interestingGroups,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates)
        .plotNovelty(
            object = metrics,
            geom = geom,
            min = min,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill,
            multiplexed = multiplexed)
    })



#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    signature("data.frame"),
    .plotNovelty)



#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    signature("seurat"),
    function(
        object,
        geom = "boxplot",
        min,
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
                .[["minNovelty"]]
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        .plotNovelty(
            object = metrics,
            geom = geom,
            min = min,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill,
            multiplexed = multiplexed)
    })
