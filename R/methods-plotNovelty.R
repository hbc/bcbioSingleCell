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
#' load(system.file("extdata/bcb_small.rda", package = "bcbioSingleCell"))
#' load(system.file("extdata/seurat_small.rda", package = "bcbioSingleCell"))
#'
#' # bcbioSingleCell ====
#' plotNovelty(bcb_small)
#'
#' # seurat ====
#' plotNovelty(seurat_small)
NULL



# Constructors =================================================================
.plotNovelty <- function(
    object,
    geom = c("violin", "boxplot", "histogram", "ridgeline"),
    interestingGroups,
    min,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    geom <- match.arg(geom)
    .plotQCMetric(
        object = object,
        metricCol = "log10GenesPerUMI",
        geom = geom,
        interestingGroups = interestingGroups,
        min = min,
        fill = fill
    )
}



# Methods ======================================================================
#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    signature("bcbioSingleCell"),
    .plotNovelty
)



#' @rdname plotNovelty
#' @export
setMethod(
    "plotNovelty",
    signature("seurat"),
    .plotNovelty
)
