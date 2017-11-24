#' Plot Mitochondrial Counts vs. Coding Counts
#'
#' @rdname plotMitoVsCoding
#' @name plotMitoVsCoding
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotUMIsVsGenes
#'
#' @examples
#' # bcbioSingleCell
#' \dontrun{
#' plotMitoVsCoding(bcb)
#' }
#'
#' # seurat
#' \dontrun{
#' plotMitoVsCoding(seurat)
#' }
#'
#' # data.frame
#' \dontrun{
#' metrics <- metrics(bcb)
#' plotMitoVsCoding(metrics)
#' }
NULL



# Constructors ====
#' @importFrom ggplot2 facet_wrap labs
#' @importFrom viridis scale_color_viridis
.plotMitoVsCoding <- function(
    object,
    interestingGroups,
    color = scale_color_viridis(discrete = TRUE)) {
    p <- .plotQCScatterplot(object, xCol = "nCoding", yCol = "nMito")

    # Label interesting groups
    if (!missing(interestingGroups)) {
        p <- p + labs(color = paste(interestingGroups, collapse = ":\n"))
    } else {
        p <- p + labs(color = NULL, fill = NULL)
    }

    # Color palette
    if (!is.null(color)) {
        p <- p + color
    }

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (!is.null(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free")
    }

    p
}



# Methods ====
#' @rdname plotMitoVsCoding
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("bcbioSingleCell"),
    function(
        object,
        interestingGroups,
        samplesOnYAxis = TRUE,
        color = scale_color_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        metrics <- metrics(
            object,
            interestingGroups = interestingGroups)
        .plotMitoVsCoding(
            object = metrics,
            interestingGroups = interestingGroups,
            color = color)
    })



#' @rdname plotMitoVsCoding
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("data.frame"),
    .plotMitoVsCoding)



#' @rdname plotMitoVsCoding
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("seurat"),
    function(
        object,
        interestingGroups,
        color = scale_color_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        .plotMitoVsCoding(
            object = metrics,
            interestingGroups = interestingGroups,
            color = color)
    })
