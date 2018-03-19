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
#' # bcbioSingleCell ====
#' plotNovelty(bcb_small)
#'
#' # seurat ====
#' plotNovelty(seurat_small)
NULL



# Constructors =================================================================
.plotNovelty <- function(
    object,
    geom = "violin",
    min,
    interestingGroups,
    samplesOnYAxis = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    metricCol <- "log10GenesPerUMI"

    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    if (missing(min)) {
        min <- metadata(object)[["filterParams"]][["minNovelty"]]
    }

    metrics <- metrics(object, interestingGroups = interestingGroups)

    p <- .plotQCGeom(
        metrics = metrics,
        geom = geom,
        metricCol = metricCol,
        min = min
    )

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
        p <- p + .medianLabels(metrics, medianCol = metricCol, digits = 2L)
    }

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(metrics))) {
        facets <- "sampleNameAggregate"
    }
    if (is.character(facets)) {
        p <- p + facet_wrap(facets = facets, scales = "free_y")
    }

    if (isTRUE(samplesOnYAxis) && geom %in% validQCGeomFlip) {
        p <- p + coord_flip()
    }

    p
}



# Methods ======================================================================
#' @rdname plotNovelty
#' @importFrom bcbioBase interestingGroups
#' @export
setMethod(
    "plotNovelty",
    signature("bcbioSingleCell"),
    .plotNovelty
)



#' @rdname plotNovelty
#' @importFrom bcbioBase interestingGroups
#' @export
setMethod(
    "plotNovelty",
    signature("seurat"),
    .plotNovelty
)
