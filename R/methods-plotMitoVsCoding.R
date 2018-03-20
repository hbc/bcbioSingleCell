#' Plot Mitochondrial Counts vs. Coding Counts
#'
#' @rdname plotMitoVsCoding
#' @name plotMitoVsCoding
#' @family Quality Control Functions
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotUMIsVsGenes
#'
#' @examples
#' # bcbioSingleCell ====
#' plotMitoVsCoding(bcb_small)
#'
#' # seurat ====
#' plotMitoVsCoding(seurat_small)
NULL



# Constructors =================================================================
#' @importFrom bcbioBase interestingGroups
#' @importFrom ggplot2 facet_wrap labs
.plotMitoVsCoding <- function(
    object,
    interestingGroups,
    color = scale_color_viridis(discrete = TRUE)
) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }

    metrics <- metrics(object, interestingGroups)

    p <- .plotQCScatterplot(metrics, xCol = "nCoding", yCol = "nMito")

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
#' @rdname plotMitoVsCoding
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("bcbioSingleCell"),
    .plotMitoVsCoding
)



#' @rdname plotMitoVsCoding
#' @export
setMethod(
    "plotMitoVsCoding",
    signature("seurat"),
    .plotMitoVsCoding
)
