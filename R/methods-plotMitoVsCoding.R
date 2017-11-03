#' Plot Mitochondrial Counts vs. Coding Counts
#'
#' @rdname plotMitoVsCoding
#' @name plotMitoVsCoding
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
#' @importFrom viridis scale_fill_viridis
.plotMitoVsCoding <- function(
    object,
    geom = "violin",
    max = Inf,
    interestingGroups,
    multiplexed = FALSE,
    samplesOnYAxis = TRUE,
    color = scale_color_viridis(discrete = TRUE)) {
    p <- .plotQCScatterplot(
            object,
            xCol = "nCoding",
            yCol = "nMito")

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
    if (isTRUE(multiplexed) & length(unique(object[["description"]])) > 1) {
        facets <- c(facets, "description")
    }
    if (!isTRUE(aggregateReplicates) &
        "sampleNameAggregate" %in% colnames(object)) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (!is.null(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free")
    }

    p
}



# Methods ====
#' @rdname plotMitoVsCoding
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("bcbioSingleCell"),
    function(
        object,
        interestingGroups,
        filterCells = TRUE,
        aggregateReplicates = TRUE,
        samplesOnYAxis = TRUE,
        color = scale_color_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- metadata(object)[["interestingGroups"]]
        }
        multiplexed <- metadata(object)[["multiplexedFASTQ"]]
        metrics <- metrics(
            object,
            interestingGroups = interestingGroups,
            filterCells = filterCells,
            aggregateReplicates = aggregateReplicates)
        .plotMitoVsCoding(
            object = metrics,
            interestingGroups = interestingGroups,
            multiplexed = multiplexed,
            color = color)
    })



#' @rdname plotMitoVsCoding
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("data.frame"),
    .plotMitoVsCoding)



#' @rdname plotMitoVsCoding
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("seurat"),
    function(
        object,
        interestingGroups,
        multiplexed = FALSE,
        color = scale_color_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- slot(object, "misc") %>%
                .[["bcbio"]] %>%
                .[["interestingGroups"]]
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        .plotMitoVsCoding(
            object = metrics,
            interestingGroups = interestingGroups,
            multiplexed = multiplexed,
            color = color)
    })
