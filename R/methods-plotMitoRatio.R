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
#' load(system.file("extdata/bcb.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' plotMitoRatio(bcb)
#'
#' # seurat
#' plotMitoRatio(seurat)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plotMitoRatio(df)
NULL



# Constructors =================================================================
#' @importFrom viridis scale_fill_viridis
.plotMitoRatio <- function(
    object,
    geom = "violin",
    max = Inf,
    interestingGroups,
    samplesOnYAxis = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)) {
    metricCol <- "mitoRatio"
    p <- .plotQCGeom(
            object,
            geom = geom,
            metricCol = metricCol,
            max = max)

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

    # Median labels
    if (geom %in% validMedianGeom) {
        p <- p + .medianLabels(object, medianCol = metricCol, digits = 2L)
    }

    # Facets
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
#' @rdname plotMitoRatio
#' @importFrom bcbioBase interestingGroups
#' @importFrom viridis scale_fill_viridis
#' @export
setMethod(
    "plotMitoRatio",
    signature("bcbioSingleCell"),
    function(
        object,
        geom = "violin",
        max,
        interestingGroups,
        samplesOnYAxis = TRUE,
        fill = viridis::scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        if (missing(max)) {
            max <- metadata(object)[["filterParams"]][["maxMitoRatio"]]
        }
        metrics <- metrics(
            object,
            interestingGroups = interestingGroups)
        .plotMitoRatio(
            object = metrics,
            geom = geom,
            max = max,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill)
    })



#' @rdname plotMitoRatio
#' @export
setMethod(
    "plotMitoRatio",
    signature("data.frame"),
    .plotMitoRatio)



#' @rdname plotMitoRatio
#' @importFrom bcbioBase interestingGroups
#' @importFrom viridis scale_fill_viridis
#' @export
setMethod(
    "plotMitoRatio",
    signature("seurat"),
    function(
        object,
        geom = "violin",
        max,
        interestingGroups,
        samplesOnYAxis = TRUE,
        fill = viridis::scale_fill_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        if (missing(max)) {
            max <- bcbio(object)[["filterParams"]][["maxMitoRatio"]]
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        .plotMitoRatio(
            object = metrics,
            geom = geom,
            max = max,
            interestingGroups = interestingGroups,
            samplesOnYAxis = samplesOnYAxis,
            fill = fill)
    })
