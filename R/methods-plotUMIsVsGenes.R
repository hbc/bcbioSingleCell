#' Plot UMI and Gene Correlation
#'
#' @rdname plotUMIsVsGenes
#' @name plotUMIsVsGenes
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
#' load(system.file(
#'     file.path("extdata", "bcb.rda"),
#'     package = "bcbioSingleCell"))
#' load(system.file(
#'     file.path("extdata", "seurat.rda"),
#'     package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell
#' plotUMIsVsGenes(bcb)
#'
#' # seurat
#' plotUMIsVsGenes(seurat)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plotUMIsPerCell(df)
NULL



# Constructors ====
#' @importFrom ggplot2 facet_wrap labs
#' @importFrom viridis scale_color_viridis
.plotUMIsVsGenes <- function(
    object,
    interestingGroups,
    color = scale_color_viridis(discrete = TRUE)) {
    p <- .plotQCScatterplot(object, xCol = "nUMI", yCol = "nGene")

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
#' @rdname plotUMIsVsGenes
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotUMIsVsGenes",
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
        .plotUMIsVsGenes(
            object = metrics,
            interestingGroups = interestingGroups,
            color = color)
    })



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("data.frame"),
    .plotUMIsVsGenes)



#' @rdname plotUMIsVsGenes
#' @importFrom viridis scale_color_viridis
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("seurat"),
    function(
        object,
        interestingGroups,
        color = scale_color_viridis(discrete = TRUE)) {
        if (missing(interestingGroups)) {
            interestingGroups <- basejump::interestingGroups(object)
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        .plotUMIsVsGenes(
            object = metrics,
            interestingGroups = interestingGroups,
            color = color)
    })
