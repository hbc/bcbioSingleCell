#' Plot UMI and Gene Correlation
#'
#' @name plotUMIsVsGenes
#' @family Quality Control Functions
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
#' # bcbioSingleCell ====
#' plotUMIsVsGenes(bcb_small)
#'
#' # seurat ====
#' plotUMIsVsGenes(pbmc_small)
#' plotUMIsVsGenes(seurat_small)
NULL



# Constructors =================================================================
#' @importFrom ggplot2 facet_wrap labs
.plotUMIsVsGenes <- function(
    object,
    interestingGroups,
    color = scale_color_viridis(discrete = TRUE)
) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }

    metrics <- metrics(object, interestingGroups)

    p <- .plotQCScatterplot(metrics, xCol = "nUMI", yCol = "nGene")

    # Label interesting groups
    if (!missing(interestingGroups)) {
        p <- p + labs(color = paste(interestingGroups, collapse = ":\n"))
    } else {
        p <- p + labs(color = NULL, fill = NULL)
    }

    # Color palette
    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(metrics))) {
        facets <- c(facets, "sampleNameAggregate")
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free")
    }

    p
}



# Methods ======================================================================
#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("bcbioSingleCell"),
    .plotUMIsVsGenes
)



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("seurat"),
    .plotUMIsVsGenes
)
