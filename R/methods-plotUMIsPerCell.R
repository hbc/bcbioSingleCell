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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' plotUMIsPerCell(bcb)
#'
#' # seurat
#' plotUMIsPerCell(seurat)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plotUMIsPerCell(df)
NULL



# Constructors =================================================================
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

    # Median labels
    if (geom %in% validMedianGeom) {
        p <- p + .medianLabels(object, medianCol = metricCol)
    }

    # Wrap aggregated samples
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
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
        samplesOnYAxis = TRUE,
        fill = scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        if (missing(min)) {
            min <- metadata(object)[["filterParams"]][["minUMIs"]]
        }
        metrics <- metrics(
            object,
            interestingGroups = interestingGroups)
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
            interestingGroups <- bcbioBase::interestingGroups(object)
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
