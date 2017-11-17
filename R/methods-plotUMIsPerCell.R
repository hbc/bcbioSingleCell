#' Plot UMIs per Cell
#'
#' Plot the universal molecular identifiers (UMIs) per cell.
#'
#' @rdname plotUMIsPerCell
#' @name plotUMIsPerCell
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#'
#' @examples
#' # bcbioSingleCell
#' \dontrun{
#' plotUMIsPerCell(bcb)
#' }
#'
#' # seurat
#' \dontrun{
#' plotUMIsPerCell(seurat)}
#'
#' # data.frame
#' \dontrun{
#' metrics <- metrics(bcb)
#' plotUMIsPerCell(metrics)
#' }
NULL



# Constructors ====
#' @importFrom viridis scale_fill_viridis
.plotUMIsPerCell <- function(
    object,
    geom = "violin",
    min = 0,
    interestingGroups,
    samplesOnYAxis = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)) {
    metricCol <- "nUMI"
    p <- .plotQCGeom(
        object,
        geom = geom,
        metricCol = metricCol,
        min = min)

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
    if (isTRUE(.checkAggregate(object))) {
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
#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("bcbioSingleCell"),
    function(
        object,
        geom = "violin",
        min,
        interestingGroups,
        filterCells = FALSE,
        samplesOnYAxis = TRUE,
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        if (missing(min)) {
            min <- metadata(object)[["filterParams"]][["minUMIs"]]
        }
        metrics <- metrics(
            object,
            interestingGroups = interestingGroups,
            filterCells = filterCells)
        .plotUMIsPerCell(
            object = metrics,
            geom = geom,
            min = min,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill)
    })



#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("data.frame"),
    .plotUMIsPerCell)



#' @rdname plotUMIsPerCell
#' @export
setMethod(
    "plotUMIsPerCell",
    signature("seurat"),
    function(
        object,
        geom = "violin",
        min,
        interestingGroups,
        samplesOnYAxis = TRUE,
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        if (missing(min)) {
            min <- bcbio(object)[["filterParams"]][["minUMIs"]]
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        .plotUMIsPerCell(
            object = metrics,
            geom = geom,
            min = min,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill)
    })
