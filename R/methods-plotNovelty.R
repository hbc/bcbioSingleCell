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
    min,
    interestingGroups,
    flip = TRUE,
    fill = scale_fill_viridis(discrete = TRUE)
) {
    if (missing(interestingGroups)) {
        interestingGroups <- bcbioBase::interestingGroups(object)
    }
    if (missing(min)) {
        min <- metadata(object)[["filterParams"]][["minNovelty"]]
    }
    .plotQCGeom(
        metrics = metrics(object, interestingGroups = interestingGroups),
        metricCol = "log10GenesPerUMI",
        geom = geom,
        min = min,
        interestingGroups = interestingGroups,
        flip = flip,
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
