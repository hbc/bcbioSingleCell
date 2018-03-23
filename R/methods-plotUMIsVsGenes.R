#' Plot UMI and Gene Correlation
#'
#' @name plotUMIsVsGenes
#' @family Quality Control Metrics
#' @author Michael Steinbaugh, Rory Kirchner
#'
#' @inherit plotGenesPerCell
#' @inheritParams general
#'
#' @examples
#' # bcbioSingleCell ====
#' plotUMIsVsGenes(bcb_small)
#'
#' # seurat ====
#' plotUMIsVsGenes(pbmc_small)
#' plotUMIsVsGenes(seurat_small)
#'
#' # data.frame
#' df <- metrics(bcb)
#' plotUMIsPerCell(df)
NULL



# Constructors =================================================================
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
    if (is(color, "ScaleDiscrete")) {
        p <- p + color
    }

    # Facets
    facets <- NULL
    if (isTRUE(.checkAggregate(object))) {
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
    function(
        object,
        interestingGroups,
        flip = TRUE,
        color = scale_color_viridis(discrete = TRUE)
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        metrics <- metrics(
            object,
            interestingGroups = interestingGroups
        )
        .plotUMIsVsGenes(
            object = metrics,
            interestingGroups = interestingGroups,
            color = color
        )
    }
)



#' @rdname plotUMIsVsGenes
#' @export
setMethod(
    "plotUMIsVsGenes",
    signature("seurat"),
    function(
        object,
        interestingGroups,
        color = scale_color_viridis(discrete = TRUE)
    ) {
        if (missing(interestingGroups)) {
            interestingGroups <- bcbioBase::interestingGroups(object)
        }
        metrics <- metrics(object, interestingGroups = interestingGroups)
        .plotUMIsVsGenes(
            object = metrics,
            interestingGroups = interestingGroups,
            color = color
        )
    }
)
